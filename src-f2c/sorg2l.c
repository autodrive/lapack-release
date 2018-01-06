#line 1 "sorg2l.f"
/* sorg2l.f -- translated by f2c (version 20100827).
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

#line 1 "sorg2l.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by s
geqlf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORG2L + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorg2l.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorg2l.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorg2l.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORG2L( M, N, K, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       REAL               A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > SORG2L generates an m by n real matrix Q with orthonormal columns, */
/* > which is defined as the last n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by SGEQLF. */
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
/* >          A is REAL array, dimension (LDA,N) */
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by SGEQLF in the last k columns of its array */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is REAL array, dimension (N) */
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

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorg2l_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l, ii;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), slarf_(char *, integer *, integer *, doublereal *, 
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

#line 148 "sorg2l.f"
    /* Parameter adjustments */
#line 148 "sorg2l.f"
    a_dim1 = *lda;
#line 148 "sorg2l.f"
    a_offset = 1 + a_dim1;
#line 148 "sorg2l.f"
    a -= a_offset;
#line 148 "sorg2l.f"
    --tau;
#line 148 "sorg2l.f"
    --work;
#line 148 "sorg2l.f"

#line 148 "sorg2l.f"
    /* Function Body */
#line 148 "sorg2l.f"
    *info = 0;
#line 149 "sorg2l.f"
    if (*m < 0) {
#line 150 "sorg2l.f"
	*info = -1;
#line 151 "sorg2l.f"
    } else if (*n < 0 || *n > *m) {
#line 152 "sorg2l.f"
	*info = -2;
#line 153 "sorg2l.f"
    } else if (*k < 0 || *k > *n) {
#line 154 "sorg2l.f"
	*info = -3;
#line 155 "sorg2l.f"
    } else if (*lda < max(1,*m)) {
#line 156 "sorg2l.f"
	*info = -5;
#line 157 "sorg2l.f"
    }
#line 158 "sorg2l.f"
    if (*info != 0) {
#line 159 "sorg2l.f"
	i__1 = -(*info);
#line 159 "sorg2l.f"
	xerbla_("SORG2L", &i__1, (ftnlen)6);
#line 160 "sorg2l.f"
	return 0;
#line 161 "sorg2l.f"
    }

/*     Quick return if possible */

#line 165 "sorg2l.f"
    if (*n <= 0) {
#line 165 "sorg2l.f"
	return 0;
#line 165 "sorg2l.f"
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

#line 170 "sorg2l.f"
    i__1 = *n - *k;
#line 170 "sorg2l.f"
    for (j = 1; j <= i__1; ++j) {
#line 171 "sorg2l.f"
	i__2 = *m;
#line 171 "sorg2l.f"
	for (l = 1; l <= i__2; ++l) {
#line 172 "sorg2l.f"
	    a[l + j * a_dim1] = 0.;
#line 173 "sorg2l.f"
/* L10: */
#line 173 "sorg2l.f"
	}
#line 174 "sorg2l.f"
	a[*m - *n + j + j * a_dim1] = 1.;
#line 175 "sorg2l.f"
/* L20: */
#line 175 "sorg2l.f"
    }

#line 177 "sorg2l.f"
    i__1 = *k;
#line 177 "sorg2l.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "sorg2l.f"
	ii = *n - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

#line 182 "sorg2l.f"
	a[*m - *n + ii + ii * a_dim1] = 1.;
#line 183 "sorg2l.f"
	i__2 = *m - *n + ii;
#line 183 "sorg2l.f"
	i__3 = ii - 1;
#line 183 "sorg2l.f"
	slarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &
		a[a_offset], lda, &work[1], (ftnlen)4);
#line 185 "sorg2l.f"
	i__2 = *m - *n + ii - 1;
#line 185 "sorg2l.f"
	d__1 = -tau[i__];
#line 185 "sorg2l.f"
	sscal_(&i__2, &d__1, &a[ii * a_dim1 + 1], &c__1);
#line 186 "sorg2l.f"
	a[*m - *n + ii + ii * a_dim1] = 1. - tau[i__];

/*        Set A(m-k+i+1:m,n-k+i) to zero */

#line 190 "sorg2l.f"
	i__2 = *m;
#line 190 "sorg2l.f"
	for (l = *m - *n + ii + 1; l <= i__2; ++l) {
#line 191 "sorg2l.f"
	    a[l + ii * a_dim1] = 0.;
#line 192 "sorg2l.f"
/* L30: */
#line 192 "sorg2l.f"
	}
#line 193 "sorg2l.f"
/* L40: */
#line 193 "sorg2l.f"
    }
#line 194 "sorg2l.f"
    return 0;

/*     End of SORG2L */

} /* sorg2l_ */

