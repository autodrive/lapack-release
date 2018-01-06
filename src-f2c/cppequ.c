#line 1 "cppequ.f"
/* cppequ.f -- translated by f2c (version 20100827).
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

#line 1 "cppequ.f"
/* > \brief \b CPPEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download CPPEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cppequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cppequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cppequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE CPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       REAL               AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               S( * ) */
/*       COMPLEX            AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > CPPEQU computes row and column scalings intended to equilibrate a */
/* > Hermitian positive definite matrix A in packed storage and reduce */
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
/* >          AP is COMPLEX array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the Hermitian matrix A, packed */
/* >          columnwise in a linear array.  The j-th column of A is stored */
/* >          in the array AP as follows: */
/* >          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j; */
/* >          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n. */
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

/* > \date December 2016 */

/* > \ingroup complexOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int cppequ_(char *uplo, integer *n, doublecomplex *ap, 
	doublereal *s, doublereal *scond, doublereal *amax, integer *info, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, jj;
    static doublereal smin;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical upper;
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

#line 160 "cppequ.f"
    /* Parameter adjustments */
#line 160 "cppequ.f"
    --s;
#line 160 "cppequ.f"
    --ap;
#line 160 "cppequ.f"

#line 160 "cppequ.f"
    /* Function Body */
#line 160 "cppequ.f"
    *info = 0;
#line 161 "cppequ.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "cppequ.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "cppequ.f"
	*info = -1;
#line 164 "cppequ.f"
    } else if (*n < 0) {
#line 165 "cppequ.f"
	*info = -2;
#line 166 "cppequ.f"
    }
#line 167 "cppequ.f"
    if (*info != 0) {
#line 168 "cppequ.f"
	i__1 = -(*info);
#line 168 "cppequ.f"
	xerbla_("CPPEQU", &i__1, (ftnlen)6);
#line 169 "cppequ.f"
	return 0;
#line 170 "cppequ.f"
    }

/*     Quick return if possible */

#line 174 "cppequ.f"
    if (*n == 0) {
#line 175 "cppequ.f"
	*scond = 1.;
#line 176 "cppequ.f"
	*amax = 0.;
#line 177 "cppequ.f"
	return 0;
#line 178 "cppequ.f"
    }

/*     Initialize SMIN and AMAX. */

#line 182 "cppequ.f"
    s[1] = ap[1].r;
#line 183 "cppequ.f"
    smin = s[1];
#line 184 "cppequ.f"
    *amax = s[1];

#line 186 "cppequ.f"
    if (upper) {

/*        UPLO = 'U':  Upper triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 191 "cppequ.f"
	jj = 1;
#line 192 "cppequ.f"
	i__1 = *n;
#line 192 "cppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 193 "cppequ.f"
	    jj += i__;
#line 194 "cppequ.f"
	    i__2 = jj;
#line 194 "cppequ.f"
	    s[i__] = ap[i__2].r;
/* Computing MIN */
#line 195 "cppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 195 "cppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 196 "cppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 196 "cppequ.f"
	    *amax = max(d__1,d__2);
#line 197 "cppequ.f"
/* L10: */
#line 197 "cppequ.f"
	}

#line 199 "cppequ.f"
    } else {

/*        UPLO = 'L':  Lower triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 204 "cppequ.f"
	jj = 1;
#line 205 "cppequ.f"
	i__1 = *n;
#line 205 "cppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 206 "cppequ.f"
	    jj = jj + *n - i__ + 2;
#line 207 "cppequ.f"
	    i__2 = jj;
#line 207 "cppequ.f"
	    s[i__] = ap[i__2].r;
/* Computing MIN */
#line 208 "cppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 208 "cppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 209 "cppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 209 "cppequ.f"
	    *amax = max(d__1,d__2);
#line 210 "cppequ.f"
/* L20: */
#line 210 "cppequ.f"
	}
#line 211 "cppequ.f"
    }

#line 213 "cppequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 217 "cppequ.f"
	i__1 = *n;
#line 217 "cppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "cppequ.f"
	    if (s[i__] <= 0.) {
#line 219 "cppequ.f"
		*info = i__;
#line 220 "cppequ.f"
		return 0;
#line 221 "cppequ.f"
	    }
#line 222 "cppequ.f"
/* L30: */
#line 222 "cppequ.f"
	}
#line 223 "cppequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 228 "cppequ.f"
	i__1 = *n;
#line 228 "cppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "cppequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 230 "cppequ.f"
/* L40: */
#line 230 "cppequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 234 "cppequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 235 "cppequ.f"
    }
#line 236 "cppequ.f"
    return 0;

/*     End of CPPEQU */

} /* cppequ_ */

