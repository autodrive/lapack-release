#line 1 "sopgtr.f"
/* sopgtr.f -- translated by f2c (version 20100827).
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

#line 1 "sopgtr.f"
/* > \brief \b SOPGTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SOPGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sopgtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sopgtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sopgtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               AP( * ), Q( LDQ, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SOPGTR generates a real orthogonal matrix Q which is defined as the */
/* > product of n-1 elementary reflectors H(i) of order n, as returned by */
/* > SSPTRD using packed storage: */
/* > */
/* > if UPLO = 'U', Q = H(n-1) . . . H(2) H(1), */
/* > */
/* > if UPLO = 'L', Q = H(1) H(2) . . . H(n-1). */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] UPLO */
/* > \verbatim */
/* >          UPLO is CHARACTER*1 */
/* >          = 'U': Upper triangular packed storage used in previous */
/* >                 call to SSPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to SSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The order of the matrix Q. N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] AP */
/* > \verbatim */
/* >          AP is REAL array, dimension (N*(N+1)/2) */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by SSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is REAL array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is REAL array, dimension (LDQ,N) */
/* >          The N-by-N orthogonal matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDQ */
/* > \verbatim */
/* >          LDQ is INTEGER */
/* >          The leading dimension of the array Q. LDQ >= max(1,N). */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N-1) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0:  successful exit */
/* >          < 0:  if INFO = -i, the i-th argument had an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date November 2011 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sopgtr_(char *uplo, integer *n, doublereal *ap, 
	doublereal *tau, doublereal *q, integer *ldq, doublereal *work, 
	integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ij;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int sorg2l_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    sorg2r_(integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), xerbla_(char *, integer *,
	     ftnlen);


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

/*     Test the input arguments */

#line 154 "sopgtr.f"
    /* Parameter adjustments */
#line 154 "sopgtr.f"
    --ap;
#line 154 "sopgtr.f"
    --tau;
#line 154 "sopgtr.f"
    q_dim1 = *ldq;
#line 154 "sopgtr.f"
    q_offset = 1 + q_dim1;
#line 154 "sopgtr.f"
    q -= q_offset;
#line 154 "sopgtr.f"
    --work;
#line 154 "sopgtr.f"

#line 154 "sopgtr.f"
    /* Function Body */
#line 154 "sopgtr.f"
    *info = 0;
#line 155 "sopgtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 156 "sopgtr.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 157 "sopgtr.f"
	*info = -1;
#line 158 "sopgtr.f"
    } else if (*n < 0) {
#line 159 "sopgtr.f"
	*info = -2;
#line 160 "sopgtr.f"
    } else if (*ldq < max(1,*n)) {
#line 161 "sopgtr.f"
	*info = -6;
#line 162 "sopgtr.f"
    }
#line 163 "sopgtr.f"
    if (*info != 0) {
#line 164 "sopgtr.f"
	i__1 = -(*info);
#line 164 "sopgtr.f"
	xerbla_("SOPGTR", &i__1, (ftnlen)6);
#line 165 "sopgtr.f"
	return 0;
#line 166 "sopgtr.f"
    }

/*     Quick return if possible */

#line 170 "sopgtr.f"
    if (*n == 0) {
#line 170 "sopgtr.f"
	return 0;
#line 170 "sopgtr.f"
    }

#line 173 "sopgtr.f"
    if (upper) {

/*        Q was determined by a call to SSPTRD with UPLO = 'U' */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the last row and column of Q equal to those of the unit */
/*        matrix */

#line 181 "sopgtr.f"
	ij = 2;
#line 182 "sopgtr.f"
	i__1 = *n - 1;
#line 182 "sopgtr.f"
	for (j = 1; j <= i__1; ++j) {
#line 183 "sopgtr.f"
	    i__2 = j - 1;
#line 183 "sopgtr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 184 "sopgtr.f"
		q[i__ + j * q_dim1] = ap[ij];
#line 185 "sopgtr.f"
		++ij;
#line 186 "sopgtr.f"
/* L10: */
#line 186 "sopgtr.f"
	    }
#line 187 "sopgtr.f"
	    ij += 2;
#line 188 "sopgtr.f"
	    q[*n + j * q_dim1] = 0.;
#line 189 "sopgtr.f"
/* L20: */
#line 189 "sopgtr.f"
	}
#line 190 "sopgtr.f"
	i__1 = *n - 1;
#line 190 "sopgtr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "sopgtr.f"
	    q[i__ + *n * q_dim1] = 0.;
#line 192 "sopgtr.f"
/* L30: */
#line 192 "sopgtr.f"
	}
#line 193 "sopgtr.f"
	q[*n + *n * q_dim1] = 1.;

/*        Generate Q(1:n-1,1:n-1) */

#line 197 "sopgtr.f"
	i__1 = *n - 1;
#line 197 "sopgtr.f"
	i__2 = *n - 1;
#line 197 "sopgtr.f"
	i__3 = *n - 1;
#line 197 "sopgtr.f"
	sorg2l_(&i__1, &i__2, &i__3, &q[q_offset], ldq, &tau[1], &work[1], &
		iinfo);

#line 199 "sopgtr.f"
    } else {

/*        Q was determined by a call to SSPTRD with UPLO = 'L'. */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the first row and column of Q equal to those of the unit */
/*        matrix */

#line 207 "sopgtr.f"
	q[q_dim1 + 1] = 1.;
#line 208 "sopgtr.f"
	i__1 = *n;
#line 208 "sopgtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 209 "sopgtr.f"
	    q[i__ + q_dim1] = 0.;
#line 210 "sopgtr.f"
/* L40: */
#line 210 "sopgtr.f"
	}
#line 211 "sopgtr.f"
	ij = 3;
#line 212 "sopgtr.f"
	i__1 = *n;
#line 212 "sopgtr.f"
	for (j = 2; j <= i__1; ++j) {
#line 213 "sopgtr.f"
	    q[j * q_dim1 + 1] = 0.;
#line 214 "sopgtr.f"
	    i__2 = *n;
#line 214 "sopgtr.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 215 "sopgtr.f"
		q[i__ + j * q_dim1] = ap[ij];
#line 216 "sopgtr.f"
		++ij;
#line 217 "sopgtr.f"
/* L50: */
#line 217 "sopgtr.f"
	    }
#line 218 "sopgtr.f"
	    ij += 2;
#line 219 "sopgtr.f"
/* L60: */
#line 219 "sopgtr.f"
	}
#line 220 "sopgtr.f"
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

#line 224 "sopgtr.f"
	    i__1 = *n - 1;
#line 224 "sopgtr.f"
	    i__2 = *n - 1;
#line 224 "sopgtr.f"
	    i__3 = *n - 1;
#line 224 "sopgtr.f"
	    sorg2r_(&i__1, &i__2, &i__3, &q[(q_dim1 << 1) + 2], ldq, &tau[1], 
		    &work[1], &iinfo);
#line 226 "sopgtr.f"
	}
#line 227 "sopgtr.f"
    }
#line 228 "sopgtr.f"
    return 0;

/*     End of SOPGTR */

} /* sopgtr_ */

