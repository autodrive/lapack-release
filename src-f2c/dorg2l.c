#line 1 "dorg2l.f"
/* dorg2l.f -- translated by f2c (version 20100827).
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

#line 1 "dorg2l.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by s
geqlf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORG2L + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO ) */

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
/* > DORG2L generates an m by n real matrix Q with orthonormal columns, */
/* > which is defined as the last n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGEQLF. */
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
/* >          The number of columns of the matrix Q. M >= N >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in] K */
/* > \verbatim */
/* >          K is INTEGER */
/* >          The number of elementary reflectors whose product defines the */
/* >          matrix Q. N >= K >= 0. */
/* > \endverbatim */
/* > */
/* > \param[in,out] A */
/* > \verbatim */
/* >          A is DOUBLE PRECISION array, dimension (LDA,N) */
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by DGEQLF in the last k columns of its array */
/* >          argument A. */
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
/* >          reflector H(i), as returned by DGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is DOUBLE PRECISION array, dimension (N) */
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
/* Subroutine */ int dorg2l_(integer *m, integer *n, integer *k, doublereal *
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

#line 148 "dorg2l.f"
    /* Parameter adjustments */
#line 148 "dorg2l.f"
    a_dim1 = *lda;
#line 148 "dorg2l.f"
    a_offset = 1 + a_dim1;
#line 148 "dorg2l.f"
    a -= a_offset;
#line 148 "dorg2l.f"
    --tau;
#line 148 "dorg2l.f"
    --work;
#line 148 "dorg2l.f"

#line 148 "dorg2l.f"
    /* Function Body */
#line 148 "dorg2l.f"
    *info = 0;
#line 149 "dorg2l.f"
    if (*m < 0) {
#line 150 "dorg2l.f"
	*info = -1;
#line 151 "dorg2l.f"
    } else if (*n < 0 || *n > *m) {
#line 152 "dorg2l.f"
	*info = -2;
#line 153 "dorg2l.f"
    } else if (*k < 0 || *k > *n) {
#line 154 "dorg2l.f"
	*info = -3;
#line 155 "dorg2l.f"
    } else if (*lda < max(1,*m)) {
#line 156 "dorg2l.f"
	*info = -5;
#line 157 "dorg2l.f"
    }
#line 158 "dorg2l.f"
    if (*info != 0) {
#line 159 "dorg2l.f"
	i__1 = -(*info);
#line 159 "dorg2l.f"
	xerbla_("DORG2L", &i__1, (ftnlen)6);
#line 160 "dorg2l.f"
	return 0;
#line 161 "dorg2l.f"
    }

/*     Quick return if possible */

#line 165 "dorg2l.f"
    if (*n <= 0) {
#line 165 "dorg2l.f"
	return 0;
#line 165 "dorg2l.f"
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

#line 170 "dorg2l.f"
    i__1 = *n - *k;
#line 170 "dorg2l.f"
    for (j = 1; j <= i__1; ++j) {
#line 171 "dorg2l.f"
	i__2 = *m;
#line 171 "dorg2l.f"
	for (l = 1; l <= i__2; ++l) {
#line 172 "dorg2l.f"
	    a[l + j * a_dim1] = 0.;
#line 173 "dorg2l.f"
/* L10: */
#line 173 "dorg2l.f"
	}
#line 174 "dorg2l.f"
	a[*m - *n + j + j * a_dim1] = 1.;
#line 175 "dorg2l.f"
/* L20: */
#line 175 "dorg2l.f"
    }

#line 177 "dorg2l.f"
    i__1 = *k;
#line 177 "dorg2l.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "dorg2l.f"
	ii = *n - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

#line 182 "dorg2l.f"
	a[*m - *n + ii + ii * a_dim1] = 1.;
#line 183 "dorg2l.f"
	i__2 = *m - *n + ii;
#line 183 "dorg2l.f"
	i__3 = ii - 1;
#line 183 "dorg2l.f"
	dlarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &
		a[a_offset], lda, &work[1], (ftnlen)4);
#line 185 "dorg2l.f"
	i__2 = *m - *n + ii - 1;
#line 185 "dorg2l.f"
	d__1 = -tau[i__];
#line 185 "dorg2l.f"
	dscal_(&i__2, &d__1, &a[ii * a_dim1 + 1], &c__1);
#line 186 "dorg2l.f"
	a[*m - *n + ii + ii * a_dim1] = 1. - tau[i__];

/*        Set A(m-k+i+1:m,n-k+i) to zero */

#line 190 "dorg2l.f"
	i__2 = *m;
#line 190 "dorg2l.f"
	for (l = *m - *n + ii + 1; l <= i__2; ++l) {
#line 191 "dorg2l.f"
	    a[l + ii * a_dim1] = 0.;
#line 192 "dorg2l.f"
/* L30: */
#line 192 "dorg2l.f"
	}
#line 193 "dorg2l.f"
/* L40: */
#line 193 "dorg2l.f"
    }
#line 194 "dorg2l.f"
    return 0;

/*     End of DORG2L */

} /* dorg2l_ */

