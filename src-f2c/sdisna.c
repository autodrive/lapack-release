#line 1 "sdisna.f"
/* sdisna.f -- translated by f2c (version 20100827).
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

#line 1 "sdisna.f"
/* > \brief \b SDISNA */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SDISNA + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sdisna.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sdisna.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sdisna.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SDISNA( JOB, M, N, D, SEP, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          JOB */
/*       INTEGER            INFO, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               D( * ), SEP( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SDISNA computes the reciprocal condition numbers for the eigenvectors */
/* > of a real symmetric or complex Hermitian matrix or for the left or */
/* > right singular vectors of a general m-by-n matrix. The reciprocal */
/* > condition number is the 'gap' between the corresponding eigenvalue or */
/* > singular value and the nearest other one. */
/* > */
/* > The bound on the error, measured by angle in radians, in the I-th */
/* > computed vector is given by */
/* > */
/* >        SLAMCH( 'E' ) * ( ANORM / SEP( I ) ) */
/* > */
/* > where ANORM = 2-norm(A) = max( abs( D(j) ) ).  SEP(I) is not allowed */
/* > to be smaller than SLAMCH( 'E' )*ANORM in order to limit the size of */
/* > the error bound. */
/* > */
/* > SDISNA may also be used to compute error bounds for eigenvectors of */
/* > the generalized symmetric definite eigenproblem. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] JOB */
/* > \verbatim */
/* >          JOB is CHARACTER*1 */
/* >          Specifies for which problem the reciprocal condition numbers */
/* >          should be computed: */
/* >          = 'E':  the eigenvectors of a symmetric/Hermitian matrix; */
/* >          = 'L':  the left singular vectors of a general matrix; */
/* >          = 'R':  the right singular vectors of a general matrix. */
/* > \endverbatim */
/* > */
/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          If JOB = 'L' or 'R', the number of columns of the matrix, */
/* >          in which case N >= 0. Ignored if JOB = 'E'. */
/* > \endverbatim */
/* > */
/* > \param[in] D */
/* > \verbatim */
/* >          D is REAL array, dimension (M) if JOB = 'E' */
/* >                              dimension (min(M,N)) if JOB = 'L' or 'R' */
/* >          The eigenvalues (if JOB = 'E') or singular values (if JOB = */
/* >          'L' or 'R') of the matrix, in either increasing or decreasing */
/* >          order. If singular values, they must be non-negative. */
/* > \endverbatim */
/* > */
/* > \param[out] SEP */
/* > \verbatim */
/* >          SEP is REAL array, dimension (M) if JOB = 'E' */
/* >                               dimension (min(M,N)) if JOB = 'L' or 'R' */
/* >          The reciprocal condition numbers of the vectors. */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit. */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup auxOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sdisna_(char *job, integer *m, integer *n, doublereal *
	d__, doublereal *sep, integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, k;
    static doublereal eps;
    static logical decr, left, incr, sing, eigen;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal anorm;
    static logical right;
    static doublereal oldgap;
    extern doublereal slamch_(char *, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal newgap, thresh;


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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input arguments */

#line 159 "sdisna.f"
    /* Parameter adjustments */
#line 159 "sdisna.f"
    --sep;
#line 159 "sdisna.f"
    --d__;
#line 159 "sdisna.f"

#line 159 "sdisna.f"
    /* Function Body */
#line 159 "sdisna.f"
    *info = 0;
#line 160 "sdisna.f"
    eigen = lsame_(job, "E", (ftnlen)1, (ftnlen)1);
#line 161 "sdisna.f"
    left = lsame_(job, "L", (ftnlen)1, (ftnlen)1);
#line 162 "sdisna.f"
    right = lsame_(job, "R", (ftnlen)1, (ftnlen)1);
#line 163 "sdisna.f"
    sing = left || right;
#line 164 "sdisna.f"
    if (eigen) {
#line 165 "sdisna.f"
	k = *m;
#line 166 "sdisna.f"
    } else if (sing) {
#line 167 "sdisna.f"
	k = min(*m,*n);
#line 168 "sdisna.f"
    }
#line 169 "sdisna.f"
    if (! eigen && ! sing) {
#line 170 "sdisna.f"
	*info = -1;
#line 171 "sdisna.f"
    } else if (*m < 0) {
#line 172 "sdisna.f"
	*info = -2;
#line 173 "sdisna.f"
    } else if (k < 0) {
#line 174 "sdisna.f"
	*info = -3;
#line 175 "sdisna.f"
    } else {
#line 176 "sdisna.f"
	incr = TRUE_;
#line 177 "sdisna.f"
	decr = TRUE_;
#line 178 "sdisna.f"
	i__1 = k - 1;
#line 178 "sdisna.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 179 "sdisna.f"
	    if (incr) {
#line 179 "sdisna.f"
		incr = incr && d__[i__] <= d__[i__ + 1];
#line 179 "sdisna.f"
	    }
#line 181 "sdisna.f"
	    if (decr) {
#line 181 "sdisna.f"
		decr = decr && d__[i__] >= d__[i__ + 1];
#line 181 "sdisna.f"
	    }
#line 183 "sdisna.f"
/* L10: */
#line 183 "sdisna.f"
	}
#line 184 "sdisna.f"
	if (sing && k > 0) {
#line 185 "sdisna.f"
	    if (incr) {
#line 185 "sdisna.f"
		incr = incr && 0. <= d__[1];
#line 185 "sdisna.f"
	    }
#line 187 "sdisna.f"
	    if (decr) {
#line 187 "sdisna.f"
		decr = decr && d__[k] >= 0.;
#line 187 "sdisna.f"
	    }
#line 189 "sdisna.f"
	}
#line 190 "sdisna.f"
	if (! (incr || decr)) {
#line 190 "sdisna.f"
	    *info = -4;
#line 190 "sdisna.f"
	}
#line 192 "sdisna.f"
    }
#line 193 "sdisna.f"
    if (*info != 0) {
#line 194 "sdisna.f"
	i__1 = -(*info);
#line 194 "sdisna.f"
	xerbla_("SDISNA", &i__1, (ftnlen)6);
#line 195 "sdisna.f"
	return 0;
#line 196 "sdisna.f"
    }

/*     Quick return if possible */

#line 200 "sdisna.f"
    if (k == 0) {
#line 200 "sdisna.f"
	return 0;
#line 200 "sdisna.f"
    }

/*     Compute reciprocal condition numbers */

#line 205 "sdisna.f"
    if (k == 1) {
#line 206 "sdisna.f"
	sep[1] = slamch_("O", (ftnlen)1);
#line 207 "sdisna.f"
    } else {
#line 208 "sdisna.f"
	oldgap = (d__1 = d__[2] - d__[1], abs(d__1));
#line 209 "sdisna.f"
	sep[1] = oldgap;
#line 210 "sdisna.f"
	i__1 = k - 1;
#line 210 "sdisna.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 211 "sdisna.f"
	    newgap = (d__1 = d__[i__ + 1] - d__[i__], abs(d__1));
#line 212 "sdisna.f"
	    sep[i__] = min(oldgap,newgap);
#line 213 "sdisna.f"
	    oldgap = newgap;
#line 214 "sdisna.f"
/* L20: */
#line 214 "sdisna.f"
	}
#line 215 "sdisna.f"
	sep[k] = oldgap;
#line 216 "sdisna.f"
    }
#line 217 "sdisna.f"
    if (sing) {
#line 218 "sdisna.f"
	if (left && *m > *n || right && *m < *n) {
#line 219 "sdisna.f"
	    if (incr) {
#line 219 "sdisna.f"
		sep[1] = min(sep[1],d__[1]);
#line 219 "sdisna.f"
	    }
#line 221 "sdisna.f"
	    if (decr) {
/* Computing MIN */
#line 221 "sdisna.f"
		d__1 = sep[k], d__2 = d__[k];
#line 221 "sdisna.f"
		sep[k] = min(d__1,d__2);
#line 221 "sdisna.f"
	    }
#line 223 "sdisna.f"
	}
#line 224 "sdisna.f"
    }

/*     Ensure that reciprocal condition numbers are not less than */
/*     threshold, in order to limit the size of the error bound */

#line 229 "sdisna.f"
    eps = slamch_("E", (ftnlen)1);
#line 230 "sdisna.f"
    safmin = slamch_("S", (ftnlen)1);
/* Computing MAX */
#line 231 "sdisna.f"
    d__2 = abs(d__[1]), d__3 = (d__1 = d__[k], abs(d__1));
#line 231 "sdisna.f"
    anorm = max(d__2,d__3);
#line 232 "sdisna.f"
    if (anorm == 0.) {
#line 233 "sdisna.f"
	thresh = eps;
#line 234 "sdisna.f"
    } else {
/* Computing MAX */
#line 235 "sdisna.f"
	d__1 = eps * anorm;
#line 235 "sdisna.f"
	thresh = max(d__1,safmin);
#line 236 "sdisna.f"
    }
#line 237 "sdisna.f"
    i__1 = k;
#line 237 "sdisna.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 238 "sdisna.f"
	d__1 = sep[i__];
#line 238 "sdisna.f"
	sep[i__] = max(d__1,thresh);
#line 239 "sdisna.f"
/* L30: */
#line 239 "sdisna.f"
    }

#line 241 "sdisna.f"
    return 0;

/*     End of SDISNA */

} /* sdisna_ */

