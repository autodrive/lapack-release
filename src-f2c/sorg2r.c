#line 1 "sorg2r.f"
/* sorg2r.f -- translated by f2c (version 20100827).
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

#line 1 "sorg2r.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b SORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by s
geqrf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download SORG2R + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorg2r.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorg2r.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorg2r.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO ) */

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
/* > SORG2R generates an m by n real matrix Q with orthonormal columns, */
/* > which is defined as the first n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* >       Q  =  H(1) H(2) . . . H(k) */
/* > */
/* > as returned by SGEQRF. */
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
/* >          On entry, the i-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by SGEQRF in the first k columns of its array */
/* >          argument A. */
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
/* >          TAU is REAL array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by SGEQRF. */
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

/* > \date September 2012 */

/* > \ingroup realOTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int sorg2r_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), slarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


/*  -- LAPACK computational routine (version 3.4.2) -- */
/*  -- LAPACK is a software package provided by Univ. of Tennessee,    -- */
/*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..-- */
/*     September 2012 */

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

#line 148 "sorg2r.f"
    /* Parameter adjustments */
#line 148 "sorg2r.f"
    a_dim1 = *lda;
#line 148 "sorg2r.f"
    a_offset = 1 + a_dim1;
#line 148 "sorg2r.f"
    a -= a_offset;
#line 148 "sorg2r.f"
    --tau;
#line 148 "sorg2r.f"
    --work;
#line 148 "sorg2r.f"

#line 148 "sorg2r.f"
    /* Function Body */
#line 148 "sorg2r.f"
    *info = 0;
#line 149 "sorg2r.f"
    if (*m < 0) {
#line 150 "sorg2r.f"
	*info = -1;
#line 151 "sorg2r.f"
    } else if (*n < 0 || *n > *m) {
#line 152 "sorg2r.f"
	*info = -2;
#line 153 "sorg2r.f"
    } else if (*k < 0 || *k > *n) {
#line 154 "sorg2r.f"
	*info = -3;
#line 155 "sorg2r.f"
    } else if (*lda < max(1,*m)) {
#line 156 "sorg2r.f"
	*info = -5;
#line 157 "sorg2r.f"
    }
#line 158 "sorg2r.f"
    if (*info != 0) {
#line 159 "sorg2r.f"
	i__1 = -(*info);
#line 159 "sorg2r.f"
	xerbla_("SORG2R", &i__1, (ftnlen)6);
#line 160 "sorg2r.f"
	return 0;
#line 161 "sorg2r.f"
    }

/*     Quick return if possible */

#line 165 "sorg2r.f"
    if (*n <= 0) {
#line 165 "sorg2r.f"
	return 0;
#line 165 "sorg2r.f"
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

#line 170 "sorg2r.f"
    i__1 = *n;
#line 170 "sorg2r.f"
    for (j = *k + 1; j <= i__1; ++j) {
#line 171 "sorg2r.f"
	i__2 = *m;
#line 171 "sorg2r.f"
	for (l = 1; l <= i__2; ++l) {
#line 172 "sorg2r.f"
	    a[l + j * a_dim1] = 0.;
#line 173 "sorg2r.f"
/* L10: */
#line 173 "sorg2r.f"
	}
#line 174 "sorg2r.f"
	a[j + j * a_dim1] = 1.;
#line 175 "sorg2r.f"
/* L20: */
#line 175 "sorg2r.f"
    }

#line 177 "sorg2r.f"
    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i) to A(i:m,i:n) from the left */

#line 181 "sorg2r.f"
	if (i__ < *n) {
#line 182 "sorg2r.f"
	    a[i__ + i__ * a_dim1] = 1.;
#line 183 "sorg2r.f"
	    i__1 = *m - i__ + 1;
#line 183 "sorg2r.f"
	    i__2 = *n - i__;
#line 183 "sorg2r.f"
	    slarf_("Left", &i__1, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[
		    i__], &a[i__ + (i__ + 1) * a_dim1], lda, &work[1], (
		    ftnlen)4);
#line 185 "sorg2r.f"
	}
#line 186 "sorg2r.f"
	if (i__ < *m) {
#line 186 "sorg2r.f"
	    i__1 = *m - i__;
#line 186 "sorg2r.f"
	    d__1 = -tau[i__];
#line 186 "sorg2r.f"
	    sscal_(&i__1, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 186 "sorg2r.f"
	}
#line 188 "sorg2r.f"
	a[i__ + i__ * a_dim1] = 1. - tau[i__];

/*        Set A(1:i-1,i) to zero */

#line 192 "sorg2r.f"
	i__1 = i__ - 1;
#line 192 "sorg2r.f"
	for (l = 1; l <= i__1; ++l) {
#line 193 "sorg2r.f"
	    a[l + i__ * a_dim1] = 0.;
#line 194 "sorg2r.f"
/* L30: */
#line 194 "sorg2r.f"
	}
#line 195 "sorg2r.f"
/* L40: */
#line 195 "sorg2r.f"
    }
#line 196 "sorg2r.f"
    return 0;

/*     End of SORG2R */

} /* sorg2r_ */

