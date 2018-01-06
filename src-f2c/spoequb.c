#line 1 "spoequb.f"
/* spoequb.f -- translated by f2c (version 20100827).
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

#line 1 "spoequb.f"
/* > \brief \b SPOEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SPOEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/spoequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/spoequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/spoequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO ) */

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
/* Subroutine */ int spoequb_(integer *n, doublereal *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal tmp, base, smin;
    extern doublereal slamch_(char *, ftnlen);
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

/*     Positive definite only performs 1 pass of equilibration. */

#line 154 "spoequb.f"
    /* Parameter adjustments */
#line 154 "spoequb.f"
    a_dim1 = *lda;
#line 154 "spoequb.f"
    a_offset = 1 + a_dim1;
#line 154 "spoequb.f"
    a -= a_offset;
#line 154 "spoequb.f"
    --s;
#line 154 "spoequb.f"

#line 154 "spoequb.f"
    /* Function Body */
#line 154 "spoequb.f"
    *info = 0;
#line 155 "spoequb.f"
    if (*n < 0) {
#line 156 "spoequb.f"
	*info = -1;
#line 157 "spoequb.f"
    } else if (*lda < max(1,*n)) {
#line 158 "spoequb.f"
	*info = -3;
#line 159 "spoequb.f"
    }
#line 160 "spoequb.f"
    if (*info != 0) {
#line 161 "spoequb.f"
	i__1 = -(*info);
#line 161 "spoequb.f"
	xerbla_("SPOEQUB", &i__1, (ftnlen)7);
#line 162 "spoequb.f"
	return 0;
#line 163 "spoequb.f"
    }

/*     Quick return if possible. */

#line 167 "spoequb.f"
    if (*n == 0) {
#line 168 "spoequb.f"
	*scond = 1.;
#line 169 "spoequb.f"
	*amax = 0.;
#line 170 "spoequb.f"
	return 0;
#line 171 "spoequb.f"
    }
#line 173 "spoequb.f"
    base = slamch_("B", (ftnlen)1);
#line 174 "spoequb.f"
    tmp = -.5 / log(base);

/*     Find the minimum and maximum diagonal elements. */

#line 178 "spoequb.f"
    s[1] = a[a_dim1 + 1];
#line 179 "spoequb.f"
    smin = s[1];
#line 180 "spoequb.f"
    *amax = s[1];
#line 181 "spoequb.f"
    i__1 = *n;
#line 181 "spoequb.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 182 "spoequb.f"
	s[i__] = a[i__ + i__ * a_dim1];
/* Computing MIN */
#line 183 "spoequb.f"
	d__1 = smin, d__2 = s[i__];
#line 183 "spoequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 184 "spoequb.f"
	d__1 = *amax, d__2 = s[i__];
#line 184 "spoequb.f"
	*amax = max(d__1,d__2);
#line 185 "spoequb.f"
/* L10: */
#line 185 "spoequb.f"
    }

#line 187 "spoequb.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 191 "spoequb.f"
	i__1 = *n;
#line 191 "spoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 192 "spoequb.f"
	    if (s[i__] <= 0.) {
#line 193 "spoequb.f"
		*info = i__;
#line 194 "spoequb.f"
		return 0;
#line 195 "spoequb.f"
	    }
#line 196 "spoequb.f"
/* L20: */
#line 196 "spoequb.f"
	}
#line 197 "spoequb.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 202 "spoequb.f"
	i__1 = *n;
#line 202 "spoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 203 "spoequb.f"
	    i__2 = (integer) (tmp * log(s[i__]));
#line 203 "spoequb.f"
	    s[i__] = pow_di(&base, &i__2);
#line 204 "spoequb.f"
/* L30: */
#line 204 "spoequb.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)). */

#line 208 "spoequb.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 209 "spoequb.f"
    }

#line 211 "spoequb.f"
    return 0;

/*     End of SPOEQUB */

} /* spoequb_ */

