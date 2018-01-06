#line 1 "dopgtr.f"
/* dopgtr.f -- translated by f2c (version 20100827).
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

#line 1 "dopgtr.f"
/* > \brief \b DOPGTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DOPGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dopgtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dopgtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dopgtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DOPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   AP( * ), Q( LDQ, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DOPGTR generates a real orthogonal matrix Q which is defined as the */
/* > product of n-1 elementary reflectors H(i) of order n, as returned by */
/* > DSPTRD using packed storage: */
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
/* >                 call to DSPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to DSPTRD. */
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
/* >          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2) */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DSPTRD. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is DOUBLE PRECISION array, dimension (LDQ,N) */
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
/* >          WORK is DOUBLE PRECISION array, dimension (N-1) */
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

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dopgtr_(char *uplo, integer *n, doublereal *ap, 
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
    extern /* Subroutine */ int dorg2l_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dorg2r_(integer *, integer *, integer *, doublereal *, integer *, 
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

#line 154 "dopgtr.f"
    /* Parameter adjustments */
#line 154 "dopgtr.f"
    --ap;
#line 154 "dopgtr.f"
    --tau;
#line 154 "dopgtr.f"
    q_dim1 = *ldq;
#line 154 "dopgtr.f"
    q_offset = 1 + q_dim1;
#line 154 "dopgtr.f"
    q -= q_offset;
#line 154 "dopgtr.f"
    --work;
#line 154 "dopgtr.f"

#line 154 "dopgtr.f"
    /* Function Body */
#line 154 "dopgtr.f"
    *info = 0;
#line 155 "dopgtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 156 "dopgtr.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 157 "dopgtr.f"
	*info = -1;
#line 158 "dopgtr.f"
    } else if (*n < 0) {
#line 159 "dopgtr.f"
	*info = -2;
#line 160 "dopgtr.f"
    } else if (*ldq < max(1,*n)) {
#line 161 "dopgtr.f"
	*info = -6;
#line 162 "dopgtr.f"
    }
#line 163 "dopgtr.f"
    if (*info != 0) {
#line 164 "dopgtr.f"
	i__1 = -(*info);
#line 164 "dopgtr.f"
	xerbla_("DOPGTR", &i__1, (ftnlen)6);
#line 165 "dopgtr.f"
	return 0;
#line 166 "dopgtr.f"
    }

/*     Quick return if possible */

#line 170 "dopgtr.f"
    if (*n == 0) {
#line 170 "dopgtr.f"
	return 0;
#line 170 "dopgtr.f"
    }

#line 173 "dopgtr.f"
    if (upper) {

/*        Q was determined by a call to DSPTRD with UPLO = 'U' */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the last row and column of Q equal to those of the unit */
/*        matrix */

#line 181 "dopgtr.f"
	ij = 2;
#line 182 "dopgtr.f"
	i__1 = *n - 1;
#line 182 "dopgtr.f"
	for (j = 1; j <= i__1; ++j) {
#line 183 "dopgtr.f"
	    i__2 = j - 1;
#line 183 "dopgtr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 184 "dopgtr.f"
		q[i__ + j * q_dim1] = ap[ij];
#line 185 "dopgtr.f"
		++ij;
#line 186 "dopgtr.f"
/* L10: */
#line 186 "dopgtr.f"
	    }
#line 187 "dopgtr.f"
	    ij += 2;
#line 188 "dopgtr.f"
	    q[*n + j * q_dim1] = 0.;
#line 189 "dopgtr.f"
/* L20: */
#line 189 "dopgtr.f"
	}
#line 190 "dopgtr.f"
	i__1 = *n - 1;
#line 190 "dopgtr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 191 "dopgtr.f"
	    q[i__ + *n * q_dim1] = 0.;
#line 192 "dopgtr.f"
/* L30: */
#line 192 "dopgtr.f"
	}
#line 193 "dopgtr.f"
	q[*n + *n * q_dim1] = 1.;

/*        Generate Q(1:n-1,1:n-1) */

#line 197 "dopgtr.f"
	i__1 = *n - 1;
#line 197 "dopgtr.f"
	i__2 = *n - 1;
#line 197 "dopgtr.f"
	i__3 = *n - 1;
#line 197 "dopgtr.f"
	dorg2l_(&i__1, &i__2, &i__3, &q[q_offset], ldq, &tau[1], &work[1], &
		iinfo);

#line 199 "dopgtr.f"
    } else {

/*        Q was determined by a call to DSPTRD with UPLO = 'L'. */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the first row and column of Q equal to those of the unit */
/*        matrix */

#line 207 "dopgtr.f"
	q[q_dim1 + 1] = 1.;
#line 208 "dopgtr.f"
	i__1 = *n;
#line 208 "dopgtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 209 "dopgtr.f"
	    q[i__ + q_dim1] = 0.;
#line 210 "dopgtr.f"
/* L40: */
#line 210 "dopgtr.f"
	}
#line 211 "dopgtr.f"
	ij = 3;
#line 212 "dopgtr.f"
	i__1 = *n;
#line 212 "dopgtr.f"
	for (j = 2; j <= i__1; ++j) {
#line 213 "dopgtr.f"
	    q[j * q_dim1 + 1] = 0.;
#line 214 "dopgtr.f"
	    i__2 = *n;
#line 214 "dopgtr.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 215 "dopgtr.f"
		q[i__ + j * q_dim1] = ap[ij];
#line 216 "dopgtr.f"
		++ij;
#line 217 "dopgtr.f"
/* L50: */
#line 217 "dopgtr.f"
	    }
#line 218 "dopgtr.f"
	    ij += 2;
#line 219 "dopgtr.f"
/* L60: */
#line 219 "dopgtr.f"
	}
#line 220 "dopgtr.f"
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

#line 224 "dopgtr.f"
	    i__1 = *n - 1;
#line 224 "dopgtr.f"
	    i__2 = *n - 1;
#line 224 "dopgtr.f"
	    i__3 = *n - 1;
#line 224 "dopgtr.f"
	    dorg2r_(&i__1, &i__2, &i__3, &q[(q_dim1 << 1) + 2], ldq, &tau[1], 
		    &work[1], &iinfo);
#line 226 "dopgtr.f"
	}
#line 227 "dopgtr.f"
    }
#line 228 "dopgtr.f"
    return 0;

/*     End of DOPGTR */

} /* dopgtr_ */

