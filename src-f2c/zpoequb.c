#line 1 "zpoequb.f"
/* zpoequb.f -- translated by f2c (version 20100827).
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

#line 1 "zpoequb.f"
/* > \brief \b ZPOEQUB */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPOEQUB + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpoequb
.f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpoequb
.f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpoequb
.f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPOEQUB( N, A, LDA, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, LDA, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ) */
/*       DOUBLE PRECISION   S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPOEQUB computes row and column scalings intended to equilibrate a */
/* > Hermitian positive definite matrix A and reduce its condition number */
/* > (with respect to the two-norm).  S contains the scale factors, */
/* > S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with */
/* > elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This */
/* > choice of S puts the condition number of B within a factor N of the */
/* > smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > */
/* > This routine differs from ZPOEQU by restricting the scaling factors */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
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

/* > \ingroup complex16POcomputational */

/*  ===================================================================== */
/* Subroutine */ int zpoequb_(integer *n, doublecomplex *a, integer *lda, 
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

#line 162 "zpoequb.f"
    /* Parameter adjustments */
#line 162 "zpoequb.f"
    a_dim1 = *lda;
#line 162 "zpoequb.f"
    a_offset = 1 + a_dim1;
#line 162 "zpoequb.f"
    a -= a_offset;
#line 162 "zpoequb.f"
    --s;
#line 162 "zpoequb.f"

#line 162 "zpoequb.f"
    /* Function Body */
#line 162 "zpoequb.f"
    *info = 0;
#line 163 "zpoequb.f"
    if (*n < 0) {
#line 164 "zpoequb.f"
	*info = -1;
#line 165 "zpoequb.f"
    } else if (*lda < max(1,*n)) {
#line 166 "zpoequb.f"
	*info = -3;
#line 167 "zpoequb.f"
    }
#line 168 "zpoequb.f"
    if (*info != 0) {
#line 169 "zpoequb.f"
	i__1 = -(*info);
#line 169 "zpoequb.f"
	xerbla_("ZPOEQUB", &i__1, (ftnlen)7);
#line 170 "zpoequb.f"
	return 0;
#line 171 "zpoequb.f"
    }

/*     Quick return if possible. */

#line 175 "zpoequb.f"
    if (*n == 0) {
#line 176 "zpoequb.f"
	*scond = 1.;
#line 177 "zpoequb.f"
	*amax = 0.;
#line 178 "zpoequb.f"
	return 0;
#line 179 "zpoequb.f"
    }
#line 181 "zpoequb.f"
    base = dlamch_("B", (ftnlen)1);
#line 182 "zpoequb.f"
    tmp = -.5 / log(base);

/*     Find the minimum and maximum diagonal elements. */

#line 186 "zpoequb.f"
    i__1 = a_dim1 + 1;
#line 186 "zpoequb.f"
    s[1] = a[i__1].r;
#line 187 "zpoequb.f"
    smin = s[1];
#line 188 "zpoequb.f"
    *amax = s[1];
#line 189 "zpoequb.f"
    i__1 = *n;
#line 189 "zpoequb.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 190 "zpoequb.f"
	i__2 = i__;
#line 190 "zpoequb.f"
	i__3 = i__ + i__ * a_dim1;
#line 190 "zpoequb.f"
	s[i__2] = a[i__3].r;
/* Computing MIN */
#line 191 "zpoequb.f"
	d__1 = smin, d__2 = s[i__];
#line 191 "zpoequb.f"
	smin = min(d__1,d__2);
/* Computing MAX */
#line 192 "zpoequb.f"
	d__1 = *amax, d__2 = s[i__];
#line 192 "zpoequb.f"
	*amax = max(d__1,d__2);
#line 193 "zpoequb.f"
/* L10: */
#line 193 "zpoequb.f"
    }

#line 195 "zpoequb.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 199 "zpoequb.f"
	i__1 = *n;
#line 199 "zpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "zpoequb.f"
	    if (s[i__] <= 0.) {
#line 201 "zpoequb.f"
		*info = i__;
#line 202 "zpoequb.f"
		return 0;
#line 203 "zpoequb.f"
	    }
#line 204 "zpoequb.f"
/* L20: */
#line 204 "zpoequb.f"
	}
#line 205 "zpoequb.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 210 "zpoequb.f"
	i__1 = *n;
#line 210 "zpoequb.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "zpoequb.f"
	    i__2 = (integer) (tmp * log(s[i__]));
#line 211 "zpoequb.f"
	    s[i__] = pow_di(&base, &i__2);
#line 212 "zpoequb.f"
/* L30: */
#line 212 "zpoequb.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)). */

#line 216 "zpoequb.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 217 "zpoequb.f"
    }

#line 219 "zpoequb.f"
    return 0;

/*     End of ZPOEQUB */

} /* zpoequb_ */

