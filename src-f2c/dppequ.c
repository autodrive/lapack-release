#line 1 "dppequ.f"
/* dppequ.f -- translated by f2c (version 20100827).
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

#line 1 "dppequ.f"
/* > \brief \b DPPEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DPPEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dppequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dppequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dppequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), S( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DPPEQU computes row and column scalings intended to equilibrate a */
/* > symmetric positive definite matrix A in packed storage and reduce */
/* > its condition number (with respect to the two-norm).  S contains the */
/* > scale factors, S(i)=1/sqrt(A(i,i)), chosen so that the scaled matrix */
/* > B with elements B(i,j)=S(i)*A(i,j)*S(j) has ones on the diagonal. */
/* > This choice of S puts the condition number of B within a factor N of */
/* > the smallest possible condition number over all possible diagonal */
/* > scalings. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U':  Upper triangle of A is stored; */
/* >          = 'L':  Lower triangle of A is stored. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix A.  N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the symmetric matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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

/* > \date November 2011 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dppequ_(char *uplo, integer *n, doublereal *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, jj;
    static doublereal smin;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
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

#line 158 "dppequ.f"
    /* Parameter adjustments */
#line 158 "dppequ.f"
    --s;
#line 158 "dppequ.f"
    --ap;
#line 158 "dppequ.f"

#line 158 "dppequ.f"
    /* Function Body */
#line 158 "dppequ.f"
    *info = 0;
#line 159 "dppequ.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 160 "dppequ.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 161 "dppequ.f"
	*info = -1;
#line 162 "dppequ.f"
    } else if (*n < 0) {
#line 163 "dppequ.f"
	*info = -2;
#line 164 "dppequ.f"
    }
#line 165 "dppequ.f"
    if (*info != 0) {
#line 166 "dppequ.f"
	i__1 = -(*info);
#line 166 "dppequ.f"
	xerbla_("DPPEQU", &i__1, (ftnlen)6);
#line 167 "dppequ.f"
	return 0;
#line 168 "dppequ.f"
    }

/*     Quick return if possible */

#line 172 "dppequ.f"
    if (*n == 0) {
#line 173 "dppequ.f"
	*scond = 1.;
#line 174 "dppequ.f"
	*amax = 0.;
#line 175 "dppequ.f"
	return 0;
#line 176 "dppequ.f"
    }

/*     Initialize SMIN and AMAX. */

#line 180 "dppequ.f"
    s[1] = ap[1];
#line 181 "dppequ.f"
    smin = s[1];
#line 182 "dppequ.f"
    *amax = s[1];

#line 184 "dppequ.f"
    if (upper) {

/*        UPLO = 'U':  Upper triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 189 "dppequ.f"
	jj = 1;
#line 190 "dppequ.f"
	i__1 = *n;
#line 190 "dppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 191 "dppequ.f"
	    jj += i__;
#line 192 "dppequ.f"
	    s[i__] = ap[jj];
/* Computing MIN */
#line 193 "dppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 193 "dppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 194 "dppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 194 "dppequ.f"
	    *amax = max(d__1,d__2);
#line 195 "dppequ.f"
/* L10: */
#line 195 "dppequ.f"
	}

#line 197 "dppequ.f"
    } else {

/*        UPLO = 'L':  Lower triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 202 "dppequ.f"
	jj = 1;
#line 203 "dppequ.f"
	i__1 = *n;
#line 203 "dppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 204 "dppequ.f"
	    jj = jj + *n - i__ + 2;
#line 205 "dppequ.f"
	    s[i__] = ap[jj];
/* Computing MIN */
#line 206 "dppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 206 "dppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 207 "dppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 207 "dppequ.f"
	    *amax = max(d__1,d__2);
#line 208 "dppequ.f"
/* L20: */
#line 208 "dppequ.f"
	}
#line 209 "dppequ.f"
    }

#line 211 "dppequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 215 "dppequ.f"
	i__1 = *n;
#line 215 "dppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 216 "dppequ.f"
	    if (s[i__] <= 0.) {
#line 217 "dppequ.f"
		*info = i__;
#line 218 "dppequ.f"
		return 0;
#line 219 "dppequ.f"
	    }
#line 220 "dppequ.f"
/* L30: */
#line 220 "dppequ.f"
	}
#line 221 "dppequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 226 "dppequ.f"
	i__1 = *n;
#line 226 "dppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 227 "dppequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 228 "dppequ.f"
/* L40: */
#line 228 "dppequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 232 "dppequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 233 "dppequ.f"
    }
#line 234 "dppequ.f"
    return 0;

/*     End of DPPEQU */

} /* dppequ_ */

