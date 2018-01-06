#line 1 "zppequ.f"
/* zppequ.f -- translated by f2c (version 20100827).
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

#line 1 "zppequ.f"
/* > \brief \b ZPPEQU */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZPPEQU + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppequ.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppequ.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppequ.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZPPEQU( UPLO, N, AP, S, SCOND, AMAX, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, N */
/*       DOUBLE PRECISION   AMAX, SCOND */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   S( * ) */
/*       COMPLEX*16         AP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZPPEQU computes row and column scalings intended to equilibrate a */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The upper or lower triangle of the Hermitian matrix A, packed */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zppequ_(char *uplo, integer *n, doublecomplex *ap, 
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

#line 160 "zppequ.f"
    /* Parameter adjustments */
#line 160 "zppequ.f"
    --s;
#line 160 "zppequ.f"
    --ap;
#line 160 "zppequ.f"

#line 160 "zppequ.f"
    /* Function Body */
#line 160 "zppequ.f"
    *info = 0;
#line 161 "zppequ.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 162 "zppequ.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 163 "zppequ.f"
	*info = -1;
#line 164 "zppequ.f"
    } else if (*n < 0) {
#line 165 "zppequ.f"
	*info = -2;
#line 166 "zppequ.f"
    }
#line 167 "zppequ.f"
    if (*info != 0) {
#line 168 "zppequ.f"
	i__1 = -(*info);
#line 168 "zppequ.f"
	xerbla_("ZPPEQU", &i__1, (ftnlen)6);
#line 169 "zppequ.f"
	return 0;
#line 170 "zppequ.f"
    }

/*     Quick return if possible */

#line 174 "zppequ.f"
    if (*n == 0) {
#line 175 "zppequ.f"
	*scond = 1.;
#line 176 "zppequ.f"
	*amax = 0.;
#line 177 "zppequ.f"
	return 0;
#line 178 "zppequ.f"
    }

/*     Initialize SMIN and AMAX. */

#line 182 "zppequ.f"
    s[1] = ap[1].r;
#line 183 "zppequ.f"
    smin = s[1];
#line 184 "zppequ.f"
    *amax = s[1];

#line 186 "zppequ.f"
    if (upper) {

/*        UPLO = 'U':  Upper triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 191 "zppequ.f"
	jj = 1;
#line 192 "zppequ.f"
	i__1 = *n;
#line 192 "zppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 193 "zppequ.f"
	    jj += i__;
#line 194 "zppequ.f"
	    i__2 = jj;
#line 194 "zppequ.f"
	    s[i__] = ap[i__2].r;
/* Computing MIN */
#line 195 "zppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 195 "zppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 196 "zppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 196 "zppequ.f"
	    *amax = max(d__1,d__2);
#line 197 "zppequ.f"
/* L10: */
#line 197 "zppequ.f"
	}

#line 199 "zppequ.f"
    } else {

/*        UPLO = 'L':  Lower triangle of A is stored. */
/*        Find the minimum and maximum diagonal elements. */

#line 204 "zppequ.f"
	jj = 1;
#line 205 "zppequ.f"
	i__1 = *n;
#line 205 "zppequ.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 206 "zppequ.f"
	    jj = jj + *n - i__ + 2;
#line 207 "zppequ.f"
	    i__2 = jj;
#line 207 "zppequ.f"
	    s[i__] = ap[i__2].r;
/* Computing MIN */
#line 208 "zppequ.f"
	    d__1 = smin, d__2 = s[i__];
#line 208 "zppequ.f"
	    smin = min(d__1,d__2);
/* Computing MAX */
#line 209 "zppequ.f"
	    d__1 = *amax, d__2 = s[i__];
#line 209 "zppequ.f"
	    *amax = max(d__1,d__2);
#line 210 "zppequ.f"
/* L20: */
#line 210 "zppequ.f"
	}
#line 211 "zppequ.f"
    }

#line 213 "zppequ.f"
    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

#line 217 "zppequ.f"
	i__1 = *n;
#line 217 "zppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 218 "zppequ.f"
	    if (s[i__] <= 0.) {
#line 219 "zppequ.f"
		*info = i__;
#line 220 "zppequ.f"
		return 0;
#line 221 "zppequ.f"
	    }
#line 222 "zppequ.f"
/* L30: */
#line 222 "zppequ.f"
	}
#line 223 "zppequ.f"
    } else {

/*        Set the scale factors to the reciprocals */
/*        of the diagonal elements. */

#line 228 "zppequ.f"
	i__1 = *n;
#line 228 "zppequ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 229 "zppequ.f"
	    s[i__] = 1. / sqrt(s[i__]);
#line 230 "zppequ.f"
/* L40: */
#line 230 "zppequ.f"
	}

/*        Compute SCOND = min(S(I)) / max(S(I)) */

#line 234 "zppequ.f"
	*scond = sqrt(smin) / sqrt(*amax);
#line 235 "zppequ.f"
    }
#line 236 "zppequ.f"
    return 0;

/*     End of ZPPEQU */

} /* zppequ_ */

