#line 1 "cpoequb.f"
/* cpoequb.f -- translated by f2c (version 20100827).
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

#line 1 "cpoequb.f"
/* > \brief \b CPOEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPOEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpoequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpoequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpoequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX            A( LDA, * ) */
/*       REAL               S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPOEQUB computes row and column scalings intended to equilibrate a */
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
/* >          A is COMPLEX array, dimension (LDA,N) */
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

/* > \ingroup complexPOcomputational */

/*  ===================================================================== */
/* Subroutine */ int cpoequb_(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
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

#line 156 "cpoequb.f"
    /* Parameter adjustments */
#line 156 "cpoequb.f"
    a_dim1 = *lda;
#line 156 "cpoequb.f"
    a_offset = 1 + a_dim1;
#line 156 "cpoequb.f"
    a -= a_offset;
#line 156 "cpoequb.f"
    --s;
#line 156 "cpoequb.f"

#line 156 "cpoequb.f"
    /* Function Body */
#line 156 "cpoequb.f"
    *info = 0;
#line 157 "cpoequb.f"
    if (*n < 0) {
#line 158 "cpoequb.f"
	*info = -1;
#line 159 "cpoequb.f"
    } else if (*lda < max(1,*n)) {
#line 160 "cpoequb.f"
	*info = -3;
#line 161 "cpoequb.f"
    }
#line 162 "cpoequb.f"
    if (*info != 0) {
#line 163 "cpoequb.f"
	i__1 = -(*info);
#line 163 "cpoequb.f"
	xerbla_("CPOEQUB", &i__1, (ftnlen)7);
#line 164 "cpoequb.f"
	return 0;
#line 165 "cpoequb.f"
    }

/*     Quick return if possible. */

#line 169 "cpoequb.f"
    if (*n == 0) {
#line 170 "cpoequb.f"
	*scond = 1.;
#line 171 "cpoequb.f"
	*amax = 0.;
#line 172 "cpoequb.f"
	return 0;
#line 173 "cpoequb.f"
    }
#line 175 "cpoequb.f"
    base = slamch_("B", (ftnlen)1);
#line 176 "cpoequb.f"
    tmp = -.5 / log(base);

/*     Find the minimum and maximum diagonal elements. */

#line 180 "cpoequb.f"
    i__1 = a_dim1 + 1;
#line 180 "cpoequb.f"
    s[1] = a[i__1].r;
#line 181 "cpoequb.f"
    smin = s[1];
#line 182 "cpoequb.f"
    *amax = s[1];
#line 183 "cpoequb.f"
    i__1 = *n;
#line 183 "cpoequb.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 184 "cpoequb.f"
	i__2 = i__;
#line 184 "cpoequb.f"
	i__3 = i__ + i__ * a_dim1;
#line 184 "cpoequb.f"
	s[i__2] = a[i__3].r;
/* Computing MIN */
#line 185 "cpoequb.f"
	d__1 = smin, d__2 = s[i__];
#line 185 "cpoequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 186 "cpoequb.f"
	d__1 = *amax, d__2 = s[i__];
#line 186 "cpoequb.f"
	*amax = max(d__1,d__2);
#line 187 "cpoequb.f"
/* L10: */
#line 187 "cpoequb.f"
    }

#line 189 "cpoequb.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 193 "cpoequb.f"
	i__1 = *n;
#line 193 "cpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "cpoequb.f"
	    if (s[i__] <= 0.) {
#line 195 "cpoequb.f"
		*info = i__;
#line 196 "cpoequb.f"
		return 0;
#line 197 "cpoequb.f"
	    }
#line 198 "cpoequb.f"
/* L20: */
#line 198 "cpoequb.f"
	}
#line 199 "cpoequb.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 204 "cpoequb.f"
	i__1 = *n;
#line 204 "cpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 205 "cpoequb.f"
	    i__2 = (integer) (tmp * log(s[i__]));
#line 205 "cpoequb.f"
	    s[i__] = pow_di(&base, &i__2);
#line 206 "cpoequb.f"
/* L30: */
#line 206 "cpoequb.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)). */

#line 210 "cpoequb.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 211 "cpoequb.f"
    }

#line 213 "cpoequb.f"
    return 0;

/*     End of CPOEQUB */

} /* cpoequb_ */

