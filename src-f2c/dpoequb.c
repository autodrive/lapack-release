#line 1 "dpoequb.f"
/* dpoequb.f -- translated by f2c (version 20100827).
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

#line 1 "dpoequb.f"
/* > \brief \b DPOEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPOEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpoequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpoequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpoequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO ) */

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
/* > DPOEQUB computes row and column scalings intended to equilibrate a */
/* > symmetric positive definite matrix A and reduce its condition number */
/* > (with respect to the two-norm).  S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > */
/* > This routine differs from DPOEQU by restricting the scaling factors */
/* > to a power of the radix.  Barring over- and underflow, scaling by */
/* > these factors introduces no additional rounding errors.  However, the */
/* > scaled diagonal entries are no longer approximately 1 but lie */
/* > between sqrt(radix) and 1/sqrt(radix). */
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
/* Subroutine */ int dpoequb_(integer *n, doublereal *a, integer *lda, 
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
    extern doublereal dlamch_(char *, ftnlen);
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
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

/*     Positive definite only performs 1 pass of equilibration. */

#line 160 "dpoequb.f"
    /* Parameter adjustments */
#line 160 "dpoequb.f"
    a_dim1 = *lda;
#line 160 "dpoequb.f"
    a_offset = 1 + a_dim1;
#line 160 "dpoequb.f"
    a -= a_offset;
#line 160 "dpoequb.f"
    --s;
#line 160 "dpoequb.f"

#line 160 "dpoequb.f"
    /* Function Body */
#line 160 "dpoequb.f"
    *info = 0;
#line 161 "dpoequb.f"
    if (*n < 0) {
#line 162 "dpoequb.f"
	*info = -1;
#line 163 "dpoequb.f"
    } else if (*lda < max(1,*n)) {
#line 164 "dpoequb.f"
	*info = -3;
#line 165 "dpoequb.f"
    }
#line 166 "dpoequb.f"
    if (*info != 0) {
#line 167 "dpoequb.f"
	i__1 = -(*info);
#line 167 "dpoequb.f"
	xerbla_("DPOEQUB", &i__1, (ftnlen)7);
#line 168 "dpoequb.f"
	return 0;
#line 169 "dpoequb.f"
    }

/*     Quick return if possible. */

#line 173 "dpoequb.f"
    if (*n == 0) {
#line 174 "dpoequb.f"
	*scond = 1.;
#line 175 "dpoequb.f"
	*amax = 0.;
#line 176 "dpoequb.f"
	return 0;
#line 177 "dpoequb.f"
    }
#line 179 "dpoequb.f"
    base = dlamch_("B", (ftnlen)1);
#line 180 "dpoequb.f"
    tmp = -.5 / log(base);

/*     Find the minimum and maximum diagonal elements. */

#line 184 "dpoequb.f"
    s[1] = a[a_dim1 + 1];
#line 185 "dpoequb.f"
    smin = s[1];
#line 186 "dpoequb.f"
    *amax = s[1];
#line 187 "dpoequb.f"
    i__1 = *n;
#line 187 "dpoequb.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 188 "dpoequb.f"
	s[i__] = a[i__ + i__ * a_dim1];
/* Computing MIN */
#line 189 "dpoequb.f"
	d__1 = smin, d__2 = s[i__];
#line 189 "dpoequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 190 "dpoequb.f"
	d__1 = *amax, d__2 = s[i__];
#line 190 "dpoequb.f"
	*amax = max(d__1,d__2);
#line 191 "dpoequb.f"
/* L10: */
#line 191 "dpoequb.f"
    }

#line 193 "dpoequb.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 197 "dpoequb.f"
	i__1 = *n;
#line 197 "dpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 198 "dpoequb.f"
	    if (s[i__] <= 0.) {
#line 199 "dpoequb.f"
		*info = i__;
#line 200 "dpoequb.f"
		return 0;
#line 201 "dpoequb.f"
	    }
#line 202 "dpoequb.f"
/* L20: */
#line 202 "dpoequb.f"
	}
#line 203 "dpoequb.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 208 "dpoequb.f"
	i__1 = *n;
#line 208 "dpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 209 "dpoequb.f"
	    i__2 = (integer) (tmp * log(s[i__]));
#line 209 "dpoequb.f"
	    s[i__] = pow_di(&base, &i__2);
#line 210 "dpoequb.f"
/* L30: */
#line 210 "dpoequb.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)). */

#line 214 "dpoequb.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 215 "dpoequb.f"
    }

#line 217 "dpoequb.f"
    return 0;

/*     End of DPOEQUB */

} /* dpoequb_ */

