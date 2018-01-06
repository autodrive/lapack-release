#line 1 "dorgr2.f"
/* dorgr2.f -- translated by f2c (version 20100827).
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

#line 1 "dorgr2.f"
/* > \brief \b DORGR2 generates all or part of the orthogonal matrix Q from an RQ factorization determined by 
sgerqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGR2( M, N, K, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > DORGR2 generates an m by n real matrix Q with orthonormal rows, */
/* > which is defined as the last m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* >       Q  =  H(1) H(2) . . . H(k) */
/* > */
/* > as returned by DGERQF. */
/* > \endverbatim */

/*  Arguments: */
/*  ========== */

/* > \param[in] M */
/* > \verbatim */
/* >          M is INTEGER */
/* >          The number of rows of the matrix Q. M >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] N */
/* > \verbatim */
/* >          N is INTEGER */
/* >          The number of columns of the matrix Q. N >= M. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines the */
/* >          matrix Q. M >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the (m-k+i)-th row must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by DGERQF in the last k rows of its array argument */
/* >          A. */
/* >          On exit, the m by n matrix Q. */
/* > \endverbatim */
/* > */
/* > \param[in] LDA */
/* > \verbatim */
/* >          LDA is INTEGER */
/* >          The first dimension of the array A. LDA >= max(1,M). */
/* > \endverbatim */
/* > */
/* > \param[in] TAU */
/* > \verbatim */
/* >          TAU is DOUBLE PRECISION array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by DGERQF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (M) */
/* > \endverbatim */
/* > */
/* > \param[out] INFO */
/* > \verbatim */
/* >          INFO is INTEGER */
/* >          = 0: successful exit */
/* >          < 0: if INFO = -i, the i-th argument has an illegal value */
/* > \endverbatim */

/*  Authors: */
/*  ======== */

/* > \author Univ. of Tennessee */
/* > \author Univ. of California Berkeley */
/* > \author Univ. of Colorado Denver */
/* > \author NAG Ltd. */

/* > \date December 2016 */

/* > \ingroup doubleOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int dorgr2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l, ii;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


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

/*     Test the input arguments */

#line 148 "dorgr2.f"
    /* Parameter adjustments */
#line 148 "dorgr2.f"
    a_dim1 = *lda;
#line 148 "dorgr2.f"
    a_offset = 1 + a_dim1;
#line 148 "dorgr2.f"
    a -= a_offset;
#line 148 "dorgr2.f"
    --tau;
#line 148 "dorgr2.f"
    --work;
#line 148 "dorgr2.f"

#line 148 "dorgr2.f"
    /* Function Body */
#line 148 "dorgr2.f"
    *info = 0;
#line 149 "dorgr2.f"
    if (*m < 0) {
#line 150 "dorgr2.f"
	*info = -1;
#line 151 "dorgr2.f"
    } else if (*n < *m) {
#line 152 "dorgr2.f"
	*info = -2;
#line 153 "dorgr2.f"
    } else if (*k < 0 || *k > *m) {
#line 154 "dorgr2.f"
	*info = -3;
#line 155 "dorgr2.f"
    } else if (*lda < max(1,*m)) {
#line 156 "dorgr2.f"
	*info = -5;
#line 157 "dorgr2.f"
    }
#line 158 "dorgr2.f"
    if (*info != 0) {
#line 159 "dorgr2.f"
	i__1 = -(*info);
#line 159 "dorgr2.f"
	xerbla_("DORGR2", &i__1, (ftnlen)6);
#line 160 "dorgr2.f"
	return 0;
#line 161 "dorgr2.f"
    }

/*     Quick return if possible */

#line 165 "dorgr2.f"
    if (*m <= 0) {
#line 165 "dorgr2.f"
	return 0;
#line 165 "dorgr2.f"
    }

#line 168 "dorgr2.f"
    if (*k < *m) {

/*        Initialise rows 1:m-k to rows of the unit matrix */

#line 172 "dorgr2.f"
	i__1 = *n;
#line 172 "dorgr2.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "dorgr2.f"
	    i__2 = *m - *k;
#line 173 "dorgr2.f"
	    for (l = 1; l <= i__2; ++l) {
#line 174 "dorgr2.f"
		a[l + j * a_dim1] = 0.;
#line 175 "dorgr2.f"
/* L10: */
#line 175 "dorgr2.f"
	    }
#line 176 "dorgr2.f"
	    if (j > *n - *m && j <= *n - *k) {
#line 176 "dorgr2.f"
		a[*m - *n + j + j * a_dim1] = 1.;
#line 176 "dorgr2.f"
	    }
#line 178 "dorgr2.f"
/* L20: */
#line 178 "dorgr2.f"
	}
#line 179 "dorgr2.f"
    }

#line 181 "dorgr2.f"
    i__1 = *k;
#line 181 "dorgr2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 182 "dorgr2.f"
	ii = *m - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the right */

#line 186 "dorgr2.f"
	a[ii + (*n - *m + ii) * a_dim1] = 1.;
#line 187 "dorgr2.f"
	i__2 = ii - 1;
#line 187 "dorgr2.f"
	i__3 = *n - *m + ii;
#line 187 "dorgr2.f"
	dlarf_("Right", &i__2, &i__3, &a[ii + a_dim1], lda, &tau[i__], &a[
		a_offset], lda, &work[1], (ftnlen)5);
#line 189 "dorgr2.f"
	i__2 = *n - *m + ii - 1;
#line 189 "dorgr2.f"
	d__1 = -tau[i__];
#line 189 "dorgr2.f"
	dscal_(&i__2, &d__1, &a[ii + a_dim1], lda);
#line 190 "dorgr2.f"
	a[ii + (*n - *m + ii) * a_dim1] = 1. - tau[i__];

/*        Set A(m-k+i,n-k+i+1:n) to zero */

#line 194 "dorgr2.f"
	i__2 = *n;
#line 194 "dorgr2.f"
	for (l = *n - *m + ii + 1; l <= i__2; ++l) {
#line 195 "dorgr2.f"
	    a[ii + l * a_dim1] = 0.;
#line 196 "dorgr2.f"
/* L30: */
#line 196 "dorgr2.f"
	}
#line 197 "dorgr2.f"
/* L40: */
#line 197 "dorgr2.f"
    }
#line 198 "dorgr2.f"
    return 0;

/*     End of DORGR2 */

} /* dorgr2_ */

