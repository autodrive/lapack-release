#line 1 "cpoequ.f"
/* cpoequ.f -- translated by f2c (version 20100827).
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

#line 1 "cpoequ.f"
/* > \brief \b CPOEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPOEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cpoequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cpoequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cpoequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPOEQU( N, A, LDA, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               S( * ) */
/*       COMPLEX            A( LDA, * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPOEQU computes row and column scalings intended to equilibrate a */
/* > Hermitian positive definite matrix A and reduce its condition number */
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
/* >          The N-by-N Hermitian positive definite matrix whose scaling */
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
/* Subroutine */ int cpoequ_(integer *n, doublecomplex *a, integer *lda, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
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

#line 150 "cpoequ.f"
    /* Parameter adjustments */
#line 150 "cpoequ.f"
    a_dim1 = *lda;
#line 150 "cpoequ.f"
    a_offset = 1 + a_dim1;
#line 150 "cpoequ.f"
    a -= a_offset;
#line 150 "cpoequ.f"
    --s;
#line 150 "cpoequ.f"

#line 150 "cpoequ.f"
    /* Function Body */
#line 150 "cpoequ.f"
    *info = 0;
#line 151 "cpoequ.f"
    if (*n < 0) {
#line 152 "cpoequ.f"
	*info = -1;
#line 153 "cpoequ.f"
    } else if (*lda < max(1,*n)) {
#line 154 "cpoequ.f"
	*info = -3;
#line 155 "cpoequ.f"
    }
#line 156 "cpoequ.f"
    if (*info != 0) {
#line 157 "cpoequ.f"
	i__1 = -(*info);
#line 157 "cpoequ.f"
	xerbla_("CPOEQU", &i__1, (ftnlen)6);
#line 158 "cpoequ.f"
	return 0;
#line 159 "cpoequ.f"
    }

/*     Quick return if possible */

#line 163 "cpoequ.f"
    if (*n == 0) {
#line 164 "cpoequ.f"
	*scond = 1.;
#line 165 "cpoequ.f"
	*amax = 0.;
#line 166 "cpoequ.f"
	return 0;
#line 167 "cpoequ.f"
    }

/*     Find the minimum and maximum diagonal elements. */

#line 171 "cpoequ.f"
    i__1 = a_dim1 + 1;
#line 171 "cpoequ.f"
    s[1] = a[i__1].r;
#line 172 "cpoequ.f"
    smin = s[1];
#line 173 "cpoequ.f"
    *amax = s[1];
#line 174 "cpoequ.f"
    i__1 = *n;
#line 174 "cpoequ.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 175 "cpoequ.f"
	i__2 = i__ + i__ * a_dim1;
#line 175 "cpoequ.f"
	s[i__] = a[i__2].r;
/* Computing MIN */
#line 176 "cpoequ.f"
	d__1 = smin, d__2 = s[i__];
#line 176 "cpoequ.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 177 "cpoequ.f"
	d__1 = *amax, d__2 = s[i__];
#line 177 "cpoequ.f"
	*amax = max(d__1,d__2);
#line 178 "cpoequ.f"
/* L10: */
#line 178 "cpoequ.f"
    }

#line 180 "cpoequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 184 "cpoequ.f"
	i__1 = *n;
#line 184 "cpoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "cpoequ.f"
	    if (s[i__] <= 0.) {
#line 186 "cpoequ.f"
		*info = i__;
#line 187 "cpoequ.f"
		return 0;
#line 188 "cpoequ.f"
	    }
#line 189 "cpoequ.f"
/* L20: */
#line 189 "cpoequ.f"
	}
#line 190 "cpoequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 195 "cpoequ.f"
	i__1 = *n;
#line 195 "cpoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 196 "cpoequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 197 "cpoequ.f"
/* L30: */
#line 197 "cpoequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 201 "cpoequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 202 "cpoequ.f"
    }
#line 203 "cpoequ.f"
    return 0;

/*     End of CPOEQU */

} /* cpoequ_ */

