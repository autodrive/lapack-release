#line 1 "dorgl2.f"
/* dorgl2.f -- translated by f2c (version 20100827).
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

#line 1 "dorgl2.f"
/* > \brief \b DORGL2 */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download DORGL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgl2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgl2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgl2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO ) */

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
/* > DORGL2 generates an m by n real matrix Q with orthonormal rows, */
/* > which is defined as the first m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by DGELQF. */
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
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by DGELQF in the first k rows of its array argument A. */
/* >          On exit, the m-by-n matrix Q. */
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
/* >          reflector H(i), as returned by DGELQF. */
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
/* Subroutine */ int dorgl2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l;
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

#line 147 "dorgl2.f"
    /* Parameter adjustments */
#line 147 "dorgl2.f"
    a_dim1 = *lda;
#line 147 "dorgl2.f"
    a_offset = 1 + a_dim1;
#line 147 "dorgl2.f"
    a -= a_offset;
#line 147 "dorgl2.f"
    --tau;
#line 147 "dorgl2.f"
    --work;
#line 147 "dorgl2.f"

#line 147 "dorgl2.f"
    /* Function Body */
#line 147 "dorgl2.f"
    *info = 0;
#line 148 "dorgl2.f"
    if (*m < 0) {
#line 149 "dorgl2.f"
	*info = -1;
#line 150 "dorgl2.f"
    } else if (*n < *m) {
#line 151 "dorgl2.f"
	*info = -2;
#line 152 "dorgl2.f"
    } else if (*k < 0 || *k > *m) {
#line 153 "dorgl2.f"
	*info = -3;
#line 154 "dorgl2.f"
    } else if (*lda < max(1,*m)) {
#line 155 "dorgl2.f"
	*info = -5;
#line 156 "dorgl2.f"
    }
#line 157 "dorgl2.f"
    if (*info != 0) {
#line 158 "dorgl2.f"
	i__1 = -(*info);
#line 158 "dorgl2.f"
	xerbla_("DORGL2", &i__1, (ftnlen)6);
#line 159 "dorgl2.f"
	return 0;
#line 160 "dorgl2.f"
    }

/*     Quick return if possible */

#line 164 "dorgl2.f"
    if (*m <= 0) {
#line 164 "dorgl2.f"
	return 0;
#line 164 "dorgl2.f"
    }

#line 167 "dorgl2.f"
    if (*k < *m) {

/*        Initialise rows k+1:m to rows of the unit matrix */

#line 171 "dorgl2.f"
	i__1 = *n;
#line 171 "dorgl2.f"
	for (j = 1; j <= i__1; ++j) {
#line 172 "dorgl2.f"
	    i__2 = *m;
#line 172 "dorgl2.f"
	    for (l = *k + 1; l <= i__2; ++l) {
#line 173 "dorgl2.f"
		a[l + j * a_dim1] = 0.;
#line 174 "dorgl2.f"
/* L10: */
#line 174 "dorgl2.f"
	    }
#line 175 "dorgl2.f"
	    if (j > *k && j <= *m) {
#line 175 "dorgl2.f"
		a[j + j * a_dim1] = 1.;
#line 175 "dorgl2.f"
	    }
#line 177 "dorgl2.f"
/* L20: */
#line 177 "dorgl2.f"
	}
#line 178 "dorgl2.f"
    }

#line 180 "dorgl2.f"
    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i) to A(i:m,i:n) from the right */

#line 184 "dorgl2.f"
	if (i__ < *n) {
#line 185 "dorgl2.f"
	    if (i__ < *m) {
#line 186 "dorgl2.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 187 "dorgl2.f"
		i__1 = *m - i__;
#line 187 "dorgl2.f"
		i__2 = *n - i__ + 1;
#line 187 "dorgl2.f"
		dlarf_("Right", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, &
			tau[i__], &a[i__ + 1 + i__ * a_dim1], lda, &work[1], (
			ftnlen)5);
#line 189 "dorgl2.f"
	    }
#line 190 "dorgl2.f"
	    i__1 = *n - i__;
#line 190 "dorgl2.f"
	    d__1 = -tau[i__];
#line 190 "dorgl2.f"
	    dscal_(&i__1, &d__1, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 191 "dorgl2.f"
	}
#line 192 "dorgl2.f"
	a[i__ + i__ * a_dim1] = 1. - tau[i__];

/*        Set A(i,1:i-1) to zero */

#line 196 "dorgl2.f"
	i__1 = i__ - 1;
#line 196 "dorgl2.f"
	for (l = 1; l <= i__1; ++l) {
#line 197 "dorgl2.f"
	    a[i__ + l * a_dim1] = 0.;
#line 198 "dorgl2.f"
/* L30: */
#line 198 "dorgl2.f"
	}
#line 199 "dorgl2.f"
/* L40: */
#line 199 "dorgl2.f"
    }
#line 200 "dorgl2.f"
    return 0;

/*     End of DORGL2 */

} /* dorgl2_ */

