#line 1 "zung2l.f"
/* zung2l.f -- translated by f2c (version 20100827).
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

#line 1 "zung2l.f"
/* Table of constant values */

static integer c__1 = 1;

/* > \brief \b ZUNG2L generates all or part of the unitary matrix Q from a QL factorization determined by cgeq
lf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNG2L + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2l.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2l.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2l.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO ) */

/*       .. Scalar Arguments .. */
/*       INTEGER            INFO, K, LDA, M, N */
/*       .. */
/*       .. Array Arguments .. */
/*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * ) */
/*       .. */


/* > \par Purpose: */
/*  ============= */
/* > */
/* > \verbatim */
/* > */
/* > ZUNG2L generates an m by n complex matrix Q with orthonormal columns, */
/* > which is defined as the last n columns of a product of k elementary */
/* > reflectors of order m */
/* > */
/* >       Q  =  H(k) . . . H(2) H(1) */
/* > */
/* > as returned by ZGEQLF. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the (n-k+i)-th column must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by ZGEQLF in the last k columns of its array */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGEQLF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (N) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zung2l_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1;

    /* Local variables */
    static integer i__, j, l, ii;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);


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

#line 149 "zung2l.f"
    /* Parameter adjustments */
#line 149 "zung2l.f"
    a_dim1 = *lda;
#line 149 "zung2l.f"
    a_offset = 1 + a_dim1;
#line 149 "zung2l.f"
    a -= a_offset;
#line 149 "zung2l.f"
    --tau;
#line 149 "zung2l.f"
    --work;
#line 149 "zung2l.f"

#line 149 "zung2l.f"
    /* Function Body */
#line 149 "zung2l.f"
    *info = 0;
#line 150 "zung2l.f"
    if (*m < 0) {
#line 151 "zung2l.f"
	*info = -1;
#line 152 "zung2l.f"
    } else if (*n < 0 || *n > *m) {
#line 153 "zung2l.f"
	*info = -2;
#line 154 "zung2l.f"
    } else if (*k < 0 || *k > *n) {
#line 155 "zung2l.f"
	*info = -3;
#line 156 "zung2l.f"
    } else if (*lda < max(1,*m)) {
#line 157 "zung2l.f"
	*info = -5;
#line 158 "zung2l.f"
    }
#line 159 "zung2l.f"
    if (*info != 0) {
#line 160 "zung2l.f"
	i__1 = -(*info);
#line 160 "zung2l.f"
	xerbla_("ZUNG2L", &i__1, (ftnlen)6);
#line 161 "zung2l.f"
	return 0;
#line 162 "zung2l.f"
    }

/*     Quick return if possible */

#line 166 "zung2l.f"
    if (*n <= 0) {
#line 166 "zung2l.f"
	return 0;
#line 166 "zung2l.f"
    }

/*     Initialise columns 1:n-k to columns of the unit matrix */

#line 171 "zung2l.f"
    i__1 = *n - *k;
#line 171 "zung2l.f"
    for (j = 1; j <= i__1; ++j) {
#line 172 "zung2l.f"
	i__2 = *m;
#line 172 "zung2l.f"
	for (l = 1; l <= i__2; ++l) {
#line 173 "zung2l.f"
	    i__3 = l + j * a_dim1;
#line 173 "zung2l.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 174 "zung2l.f"
/* L10: */
#line 174 "zung2l.f"
	}
#line 175 "zung2l.f"
	i__2 = *m - *n + j + j * a_dim1;
#line 175 "zung2l.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 176 "zung2l.f"
/* L20: */
#line 176 "zung2l.f"
    }

#line 178 "zung2l.f"
    i__1 = *k;
#line 178 "zung2l.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 179 "zung2l.f"
	ii = *n - *k + i__;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left */

#line 183 "zung2l.f"
	i__2 = *m - *n + ii + ii * a_dim1;
#line 183 "zung2l.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 184 "zung2l.f"
	i__2 = *m - *n + ii;
#line 184 "zung2l.f"
	i__3 = ii - 1;
#line 184 "zung2l.f"
	zlarf_("Left", &i__2, &i__3, &a[ii * a_dim1 + 1], &c__1, &tau[i__], &
		a[a_offset], lda, &work[1], (ftnlen)4);
#line 186 "zung2l.f"
	i__2 = *m - *n + ii - 1;
#line 186 "zung2l.f"
	i__3 = i__;
#line 186 "zung2l.f"
	z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 186 "zung2l.f"
	zscal_(&i__2, &z__1, &a[ii * a_dim1 + 1], &c__1);
#line 187 "zung2l.f"
	i__2 = *m - *n + ii + ii * a_dim1;
#line 187 "zung2l.f"
	i__3 = i__;
#line 187 "zung2l.f"
	z__1.r = 1. - tau[i__3].r, z__1.i = 0. - tau[i__3].i;
#line 187 "zung2l.f"
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;

/*        Set A(m-k+i+1:m,n-k+i) to zero */

#line 191 "zung2l.f"
	i__2 = *m;
#line 191 "zung2l.f"
	for (l = *m - *n + ii + 1; l <= i__2; ++l) {
#line 192 "zung2l.f"
	    i__3 = l + ii * a_dim1;
#line 192 "zung2l.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 193 "zung2l.f"
/* L30: */
#line 193 "zung2l.f"
	}
#line 194 "zung2l.f"
/* L40: */
#line 194 "zung2l.f"
    }
#line 195 "zung2l.f"
    return 0;

/*     End of ZUNG2L */

} /* zung2l_ */

