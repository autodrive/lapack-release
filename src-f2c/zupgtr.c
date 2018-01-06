#line 1 "zupgtr.f"
/* zupgtr.f -- translated by f2c (version 20100827).
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

#line 1 "zupgtr.f"
/* > \brief \b ZUPGTR */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUPGTR + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zupgtr.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zupgtr.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zupgtr.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUPGTR( UPLO, N, AP, TAU, Q, LDQ, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       CHARACTER          UPLO */
/*       INTEGER            INFO, LDQ, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         AP( * ), Q( LDQ, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUPGTR generates a complex unitary matrix Q which is defined as the */
/* > product of n-1 elementary reflectors H(i) of order n, as returned by */
/* > ZHPTRD using packed storage: */
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
/* >                 call to ZHPTRD; */
/* >          = 'L': Lower triangular packed storage used in previous */
/* >                 call to ZHPTRD. */
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
/* >          AP is COMPLEX*16 array, dimension (N*(N+1)/2) */
/* >          The vectors which define the elementary reflectors, as */
/* >          returned by ZHPTRD. */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is COMPLEX*16 array, dimension (N-1) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZHPTRD. */
/* > \endverbatim */
/* > */
/* > \param[out] Q */
/* > \verbatim */
/* >          Q is COMPLEX*16 array, dimension (LDQ,N) */
/* >          The N-by-N unitary matrix Q. */
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
/* >          WORK is COMPLEX*16 array, dimension (N-1) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zupgtr_(char *uplo, integer *n, doublecomplex *ap, 
	doublecomplex *tau, doublecomplex *q, integer *ldq, doublecomplex *
	work, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer q_dim1, q_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, ij;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer iinfo;
    static logical upper;
    extern /* Subroutine */ int zung2l_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), zung2r_(integer *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *), xerbla_(char *, integer *, ftnlen);


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

#line 155 "zupgtr.f"
    /* Parameter adjustments */
#line 155 "zupgtr.f"
    --ap;
#line 155 "zupgtr.f"
    --tau;
#line 155 "zupgtr.f"
    q_dim1 = *ldq;
#line 155 "zupgtr.f"
    q_offset = 1 + q_dim1;
#line 155 "zupgtr.f"
    q -= q_offset;
#line 155 "zupgtr.f"
    --work;
#line 155 "zupgtr.f"

#line 155 "zupgtr.f"
    /* Function Body */
#line 155 "zupgtr.f"
    *info = 0;
#line 156 "zupgtr.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 157 "zupgtr.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 158 "zupgtr.f"
	*info = -1;
#line 159 "zupgtr.f"
    } else if (*n < 0) {
#line 160 "zupgtr.f"
	*info = -2;
#line 161 "zupgtr.f"
    } else if (*ldq < max(1,*n)) {
#line 162 "zupgtr.f"
	*info = -6;
#line 163 "zupgtr.f"
    }
#line 164 "zupgtr.f"
    if (*info != 0) {
#line 165 "zupgtr.f"
	i__1 = -(*info);
#line 165 "zupgtr.f"
	xerbla_("ZUPGTR", &i__1, (ftnlen)6);
#line 166 "zupgtr.f"
	return 0;
#line 167 "zupgtr.f"
    }

/*     Quick return if possible */

#line 171 "zupgtr.f"
    if (*n == 0) {
#line 171 "zupgtr.f"
	return 0;
#line 171 "zupgtr.f"
    }

#line 174 "zupgtr.f"
    if (upper) {

/*        Q was determined by a call to ZHPTRD with UPLO = 'U' */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the last row and column of Q equal to those of the unit */
/*        matrix */

#line 182 "zupgtr.f"
	ij = 2;
#line 183 "zupgtr.f"
	i__1 = *n - 1;
#line 183 "zupgtr.f"
	for (j = 1; j <= i__1; ++j) {
#line 184 "zupgtr.f"
	    i__2 = j - 1;
#line 184 "zupgtr.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 185 "zupgtr.f"
		i__3 = i__ + j * q_dim1;
#line 185 "zupgtr.f"
		i__4 = ij;
#line 185 "zupgtr.f"
		q[i__3].r = ap[i__4].r, q[i__3].i = ap[i__4].i;
#line 186 "zupgtr.f"
		++ij;
#line 187 "zupgtr.f"
/* L10: */
#line 187 "zupgtr.f"
	    }
#line 188 "zupgtr.f"
	    ij += 2;
#line 189 "zupgtr.f"
	    i__2 = *n + j * q_dim1;
#line 189 "zupgtr.f"
	    q[i__2].r = 0., q[i__2].i = 0.;
#line 190 "zupgtr.f"
/* L20: */
#line 190 "zupgtr.f"
	}
#line 191 "zupgtr.f"
	i__1 = *n - 1;
#line 191 "zupgtr.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 192 "zupgtr.f"
	    i__2 = i__ + *n * q_dim1;
#line 192 "zupgtr.f"
	    q[i__2].r = 0., q[i__2].i = 0.;
#line 193 "zupgtr.f"
/* L30: */
#line 193 "zupgtr.f"
	}
#line 194 "zupgtr.f"
	i__1 = *n + *n * q_dim1;
#line 194 "zupgtr.f"
	q[i__1].r = 1., q[i__1].i = 0.;

/*        Generate Q(1:n-1,1:n-1) */

#line 198 "zupgtr.f"
	i__1 = *n - 1;
#line 198 "zupgtr.f"
	i__2 = *n - 1;
#line 198 "zupgtr.f"
	i__3 = *n - 1;
#line 198 "zupgtr.f"
	zung2l_(&i__1, &i__2, &i__3, &q[q_offset], ldq, &tau[1], &work[1], &
		iinfo);

#line 200 "zupgtr.f"
    } else {

/*        Q was determined by a call to ZHPTRD with UPLO = 'L'. */

/*        Unpack the vectors which define the elementary reflectors and */
/*        set the first row and column of Q equal to those of the unit */
/*        matrix */

#line 208 "zupgtr.f"
	i__1 = q_dim1 + 1;
#line 208 "zupgtr.f"
	q[i__1].r = 1., q[i__1].i = 0.;
#line 209 "zupgtr.f"
	i__1 = *n;
#line 209 "zupgtr.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 210 "zupgtr.f"
	    i__2 = i__ + q_dim1;
#line 210 "zupgtr.f"
	    q[i__2].r = 0., q[i__2].i = 0.;
#line 211 "zupgtr.f"
/* L40: */
#line 211 "zupgtr.f"
	}
#line 212 "zupgtr.f"
	ij = 3;
#line 213 "zupgtr.f"
	i__1 = *n;
#line 213 "zupgtr.f"
	for (j = 2; j <= i__1; ++j) {
#line 214 "zupgtr.f"
	    i__2 = j * q_dim1 + 1;
#line 214 "zupgtr.f"
	    q[i__2].r = 0., q[i__2].i = 0.;
#line 215 "zupgtr.f"
	    i__2 = *n;
#line 215 "zupgtr.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 216 "zupgtr.f"
		i__3 = i__ + j * q_dim1;
#line 216 "zupgtr.f"
		i__4 = ij;
#line 216 "zupgtr.f"
		q[i__3].r = ap[i__4].r, q[i__3].i = ap[i__4].i;
#line 217 "zupgtr.f"
		++ij;
#line 218 "zupgtr.f"
/* L50: */
#line 218 "zupgtr.f"
	    }
#line 219 "zupgtr.f"
	    ij += 2;
#line 220 "zupgtr.f"
/* L60: */
#line 220 "zupgtr.f"
	}
#line 221 "zupgtr.f"
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

#line 225 "zupgtr.f"
	    i__1 = *n - 1;
#line 225 "zupgtr.f"
	    i__2 = *n - 1;
#line 225 "zupgtr.f"
	    i__3 = *n - 1;
#line 225 "zupgtr.f"
	    zung2r_(&i__1, &i__2, &i__3, &q[(q_dim1 << 1) + 2], ldq, &tau[1], 
		    &work[1], &iinfo);
#line 227 "zupgtr.f"
	}
#line 228 "zupgtr.f"
    }
#line 229 "zupgtr.f"
    return 0;

/*     End of ZUPGTR */

} /* zupgtr_ */

