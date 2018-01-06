#line 1 "dpoequ.f"
/* dpoequ.f -- translated by f2c (version 20100827).
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

#line 1 "dpoequ.f"
/* > \brief \b DPOEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpoequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpoequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpoequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOEQU( N, A, LDA, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPOEQU computes row and column scalings intended to equilibrate a */
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
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
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
/* >          S is DOUBLE PRECISION array, dimension (N) */
/* >          If INFO = 0, S contains the scale factors for A. */
/* > \endverbatim */
/* > */
/* > \param[out] SCOND */
/* > \verbatim */
/* >          SCOND is DOUBLE PRECISION */
/* >          If INFO = 0, S contains the ratio of the smallest S(i) to */
/* >          the largest S(i).  If SCOND >= 0.1 and AMAX is neither too */
/* >          large nor too small, it is not worth scaling by S. */
/* > \endverbatim */
/* > */
/* > \param[out] AMAX */
/* > \verbatim */
/* >          AMAX is DOUBLE PRECISION */
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

/* > \date December 2016 */

/* > \ingroup doublePOcomputational */

/*  ===================================================================== */
/* Subroutine */ int dpoequ_(integer *n, doublereal *a, integer *lda, 
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


/*  -- LAPACK computational routine (version 3.7.0) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     December 2016 */

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

#line 148 "dpoequ.f"
    /* Parameter adjustments */
#line 148 "dpoequ.f"
    a_dim1 = *lda;
#line 148 "dpoequ.f"
    a_offset = 1 + a_dim1;
#line 148 "dpoequ.f"
    a -= a_offset;
#line 148 "dpoequ.f"
    --s;
#line 148 "dpoequ.f"

#line 148 "dpoequ.f"
    /* Function Body */
#line 148 "dpoequ.f"
    *info = 0;
#line 149 "dpoequ.f"
    if (*n < 0) {
#line 150 "dpoequ.f"
	*info = -1;
#line 151 "dpoequ.f"
    } else if (*lda < max(1,*n)) {
#line 152 "dpoequ.f"
	*info = -3;
#line 153 "dpoequ.f"
    }
#line 154 "dpoequ.f"
    if (*info != 0) {
#line 155 "dpoequ.f"
	i__1 = -(*info);
#line 155 "dpoequ.f"
	xerbla_("DPOEQU", &i__1, (ftnlen)6);
#line 156 "dpoequ.f"
	return 0;
#line 157 "dpoequ.f"
    }

/*     Quick return if possible */

#line 161 "dpoequ.f"
    if (*n == 0) {
#line 162 "dpoequ.f"
	*scond = 1.;
#line 163 "dpoequ.f"
	*amax = 0.;
#line 164 "dpoequ.f"
	return 0;
#line 165 "dpoequ.f"
    }

/*     Find the minimum and maximum diagonal elements. */

#line 169 "dpoequ.f"
    s[1] = a[a_dim1 + 1];
#line 170 "dpoequ.f"
    smin = s[1];
#line 171 "dpoequ.f"
    *amax = s[1];
#line 172 "dpoequ.f"
    i__1 = *n;
#line 172 "dpoequ.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 173 "dpoequ.f"
	s[i__] = a[i__ + i__ * a_dim1];
/* Computing MIN */
#line 174 "dpoequ.f"
	d__1 = smin, d__2 = s[i__];
#line 174 "dpoequ.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 175 "dpoequ.f"
	d__1 = *amax, d__2 = s[i__];
#line 175 "dpoequ.f"
	*amax = max(d__1,d__2);
#line 176 "dpoequ.f"
/* L10: */
#line 176 "dpoequ.f"
    }

#line 178 "dpoequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 182 "dpoequ.f"
	i__1 = *n;
#line 182 "dpoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 183 "dpoequ.f"
	    if (s[i__] <= 0.) {
#line 184 "dpoequ.f"
		*info = i__;
#line 185 "dpoequ.f"
		return 0;
#line 186 "dpoequ.f"
	    }
#line 187 "dpoequ.f"
/* L20: */
#line 187 "dpoequ.f"
	}
#line 188 "dpoequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 193 "dpoequ.f"
	i__1 = *n;
#line 193 "dpoequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "dpoequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 195 "dpoequ.f"
/* L30: */
#line 195 "dpoequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 199 "dpoequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 200 "dpoequ.f"
    }
#line 201 "dpoequ.f"
    return 0;

/*     End of DPOEQU */

} /* dpoequ_ */

