#line 1 "zungl2.f"
/* zungl2.f -- translated by f2c (version 20100827).
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

#line 1 "zungl2.f"
/* > \brief \b ZUNGL2 generates all or part of the unitary matrix Q from an LQ factorization determined by cge
lqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGL2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungl2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGL2( M, N, K, A, LDA, TAU, WORK, INFO ) */

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
/* > ZUNGL2 generates an m-by-n complex matrix Q with orthonormal rows, */
/* > which is defined as the first m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* >       Q  =  H(k)**H . . . H(2)**H H(1)**H */
/* > */
/* > as returned by ZGELQF. */
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
/* >          A is COMPLEX*16 array, dimension (LDA,N) */
/* >          On entry, the i-th row must contain the vector which defines */
/* >          the elementary reflector H(i), for i = 1,2,...,k, as returned */
/* >          by ZGELQF in the first k rows of its array argument A. */
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
/* >          TAU is COMPLEX*16 array, dimension (K) */
/* >          TAU(i) must contain the scalar factor of the elementary */
/* >          reflector H(i), as returned by ZGELQF. */
/* > \endverbatim */
/* > */
/* > \param[out] WORK */
/* > \verbatim */
/* >          WORK is COMPLEX*16 array, dimension (M) */
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

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungl2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), zlacgv_(integer *, doublecomplex *, integer *);


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

#line 148 "zungl2.f"
    /* Parameter adjustments */
#line 148 "zungl2.f"
    a_dim1 = *lda;
#line 148 "zungl2.f"
    a_offset = 1 + a_dim1;
#line 148 "zungl2.f"
    a -= a_offset;
#line 148 "zungl2.f"
    --tau;
#line 148 "zungl2.f"
    --work;
#line 148 "zungl2.f"

#line 148 "zungl2.f"
    /* Function Body */
#line 148 "zungl2.f"
    *info = 0;
#line 149 "zungl2.f"
    if (*m < 0) {
#line 150 "zungl2.f"
	*info = -1;
#line 151 "zungl2.f"
    } else if (*n < *m) {
#line 152 "zungl2.f"
	*info = -2;
#line 153 "zungl2.f"
    } else if (*k < 0 || *k > *m) {
#line 154 "zungl2.f"
	*info = -3;
#line 155 "zungl2.f"
    } else if (*lda < max(1,*m)) {
#line 156 "zungl2.f"
	*info = -5;
#line 157 "zungl2.f"
    }
#line 158 "zungl2.f"
    if (*info != 0) {
#line 159 "zungl2.f"
	i__1 = -(*info);
#line 159 "zungl2.f"
	xerbla_("ZUNGL2", &i__1, (ftnlen)6);
#line 160 "zungl2.f"
	return 0;
#line 161 "zungl2.f"
    }

/*     Quick return if possible */

#line 165 "zungl2.f"
    if (*m <= 0) {
#line 165 "zungl2.f"
	return 0;
#line 165 "zungl2.f"
    }

#line 168 "zungl2.f"
    if (*k < *m) {

/*        Initialise rows k+1:m to rows of the unit matrix */

#line 172 "zungl2.f"
	i__1 = *n;
#line 172 "zungl2.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "zungl2.f"
	    i__2 = *m;
#line 173 "zungl2.f"
	    for (l = *k + 1; l <= i__2; ++l) {
#line 174 "zungl2.f"
		i__3 = l + j * a_dim1;
#line 174 "zungl2.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 175 "zungl2.f"
/* L10: */
#line 175 "zungl2.f"
	    }
#line 176 "zungl2.f"
	    if (j > *k && j <= *m) {
#line 176 "zungl2.f"
		i__2 = j + j * a_dim1;
#line 176 "zungl2.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 176 "zungl2.f"
	    }
#line 178 "zungl2.f"
/* L20: */
#line 178 "zungl2.f"
	}
#line 179 "zungl2.f"
    }

#line 181 "zungl2.f"
    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i)**H to A(i:m,i:n) from the right */

#line 185 "zungl2.f"
	if (i__ < *n) {
#line 186 "zungl2.f"
	    i__1 = *n - i__;
#line 186 "zungl2.f"
	    zlacgv_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 187 "zungl2.f"
	    if (i__ < *m) {
#line 188 "zungl2.f"
		i__1 = i__ + i__ * a_dim1;
#line 188 "zungl2.f"
		a[i__1].r = 1., a[i__1].i = 0.;
#line 189 "zungl2.f"
		i__1 = *m - i__;
#line 189 "zungl2.f"
		i__2 = *n - i__ + 1;
#line 189 "zungl2.f"
		d_cnjg(&z__1, &tau[i__]);
#line 189 "zungl2.f"
		zlarf_("Right", &i__1, &i__2, &a[i__ + i__ * a_dim1], lda, &
			z__1, &a[i__ + 1 + i__ * a_dim1], lda, &work[1], (
			ftnlen)5);
#line 191 "zungl2.f"
	    }
#line 192 "zungl2.f"
	    i__1 = *n - i__;
#line 192 "zungl2.f"
	    i__2 = i__;
#line 192 "zungl2.f"
	    z__1.r = -tau[i__2].r, z__1.i = -tau[i__2].i;
#line 192 "zungl2.f"
	    zscal_(&i__1, &z__1, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 193 "zungl2.f"
	    i__1 = *n - i__;
#line 193 "zungl2.f"
	    zlacgv_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 194 "zungl2.f"
	}
#line 195 "zungl2.f"
	i__1 = i__ + i__ * a_dim1;
#line 195 "zungl2.f"
	d_cnjg(&z__2, &tau[i__]);
#line 195 "zungl2.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 195 "zungl2.f"
	a[i__1].r = z__1.r, a[i__1].i = z__1.i;

/*        Set A(i,1:i-1) to zero */

#line 199 "zungl2.f"
	i__1 = i__ - 1;
#line 199 "zungl2.f"
	for (l = 1; l <= i__1; ++l) {
#line 200 "zungl2.f"
	    i__2 = i__ + l * a_dim1;
#line 200 "zungl2.f"
	    a[i__2].r = 0., a[i__2].i = 0.;
#line 201 "zungl2.f"
/* L30: */
#line 201 "zungl2.f"
	}
#line 202 "zungl2.f"
/* L40: */
#line 202 "zungl2.f"
    }
#line 203 "zungl2.f"
    return 0;

/*     End of ZUNGL2 */

} /* zungl2_ */

