#line 1 "zungr2.f"
/* zungr2.f -- translated by f2c (version 20100827).
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

#line 1 "zungr2.f"
/* > \brief \b ZUNGR2 generates all or part of the unitary matrix Q from an RQ factorization determined by cge
rqf (unblocked algorithm). */

/*  =========== DOCUMENTATION =========== */

/* Online html documentation available at */
/*            http://www.netlib.org/lapack/explore-html/ */

/* > \htmlonly */
/* > Download ZUNGR2 + dependencies */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungr2.
f"> */
/* > [TGZ]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungr2.
f"> */
/* > [ZIP]</a> */
/* > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungr2.
f"> */
/* > [TXT]</a> */
/* > \endhtmlonly */

/*  Definition: */
/*  =========== */

/*       SUBROUTINE ZUNGR2( M, N, K, A, LDA, TAU, WORK, INFO ) */

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
/* > ZUNGR2 generates an m by n complex matrix Q with orthonormal rows, */
/* > which is defined as the last m rows of a product of k elementary */
/* > reflectors of order n */
/* > */
/* >       Q  =  H(1)**H H(2)**H . . . H(k)**H */
/* > */
/* > as returned by ZGERQF. */
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
/* >          On entry, the (m-k+i)-th row must contain the vector which */
/* >          defines the elementary reflector H(i), for i = 1,2,...,k, as */
/* >          returned by ZGERQF in the last k rows of its array argument */
/* >          A. */
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
/* >          reflector H(i), as returned by ZGERQF. */
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

/* > \date December 2016 */

/* > \ingroup complex16OTHERcomputational */

/*  ===================================================================== */
/* Subroutine */ int zungr2_(integer *m, integer *n, integer *k, 
	doublecomplex *a, integer *lda, doublecomplex *tau, doublecomplex *
	work, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, l, ii;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), zlacgv_(integer *, doublecomplex *, integer *);


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

#line 149 "zungr2.f"
    /* Parameter adjustments */
#line 149 "zungr2.f"
    a_dim1 = *lda;
#line 149 "zungr2.f"
    a_offset = 1 + a_dim1;
#line 149 "zungr2.f"
    a -= a_offset;
#line 149 "zungr2.f"
    --tau;
#line 149 "zungr2.f"
    --work;
#line 149 "zungr2.f"

#line 149 "zungr2.f"
    /* Function Body */
#line 149 "zungr2.f"
    *info = 0;
#line 150 "zungr2.f"
    if (*m < 0) {
#line 151 "zungr2.f"
	*info = -1;
#line 152 "zungr2.f"
    } else if (*n < *m) {
#line 153 "zungr2.f"
	*info = -2;
#line 154 "zungr2.f"
    } else if (*k < 0 || *k > *m) {
#line 155 "zungr2.f"
	*info = -3;
#line 156 "zungr2.f"
    } else if (*lda < max(1,*m)) {
#line 157 "zungr2.f"
	*info = -5;
#line 158 "zungr2.f"
    }
#line 159 "zungr2.f"
    if (*info != 0) {
#line 160 "zungr2.f"
	i__1 = -(*info);
#line 160 "zungr2.f"
	xerbla_("ZUNGR2", &i__1, (ftnlen)6);
#line 161 "zungr2.f"
	return 0;
#line 162 "zungr2.f"
    }

/*     Quick return if possible */

#line 166 "zungr2.f"
    if (*m <= 0) {
#line 166 "zungr2.f"
	return 0;
#line 166 "zungr2.f"
    }

#line 169 "zungr2.f"
    if (*k < *m) {

/*        Initialise rows 1:m-k to rows of the unit matrix */

#line 173 "zungr2.f"
	i__1 = *n;
#line 173 "zungr2.f"
	for (j = 1; j <= i__1; ++j) {
#line 174 "zungr2.f"
	    i__2 = *m - *k;
#line 174 "zungr2.f"
	    for (l = 1; l <= i__2; ++l) {
#line 175 "zungr2.f"
		i__3 = l + j * a_dim1;
#line 175 "zungr2.f"
		a[i__3].r = 0., a[i__3].i = 0.;
#line 176 "zungr2.f"
/* L10: */
#line 176 "zungr2.f"
	    }
#line 177 "zungr2.f"
	    if (j > *n - *m && j <= *n - *k) {
#line 177 "zungr2.f"
		i__2 = *m - *n + j + j * a_dim1;
#line 177 "zungr2.f"
		a[i__2].r = 1., a[i__2].i = 0.;
#line 177 "zungr2.f"
	    }
#line 179 "zungr2.f"
/* L20: */
#line 179 "zungr2.f"
	}
#line 180 "zungr2.f"
    }

#line 182 "zungr2.f"
    i__1 = *k;
#line 182 "zungr2.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 183 "zungr2.f"
	ii = *m - *k + i__;

/*        Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right */

#line 187 "zungr2.f"
	i__2 = *n - *m + ii - 1;
#line 187 "zungr2.f"
	zlacgv_(&i__2, &a[ii + a_dim1], lda);
#line 188 "zungr2.f"
	i__2 = ii + (*n - *m + ii) * a_dim1;
#line 188 "zungr2.f"
	a[i__2].r = 1., a[i__2].i = 0.;
#line 189 "zungr2.f"
	i__2 = ii - 1;
#line 189 "zungr2.f"
	i__3 = *n - *m + ii;
#line 189 "zungr2.f"
	d_cnjg(&z__1, &tau[i__]);
#line 189 "zungr2.f"
	zlarf_("Right", &i__2, &i__3, &a[ii + a_dim1], lda, &z__1, &a[
		a_offset], lda, &work[1], (ftnlen)5);
#line 191 "zungr2.f"
	i__2 = *n - *m + ii - 1;
#line 191 "zungr2.f"
	i__3 = i__;
#line 191 "zungr2.f"
	z__1.r = -tau[i__3].r, z__1.i = -tau[i__3].i;
#line 191 "zungr2.f"
	zscal_(&i__2, &z__1, &a[ii + a_dim1], lda);
#line 192 "zungr2.f"
	i__2 = *n - *m + ii - 1;
#line 192 "zungr2.f"
	zlacgv_(&i__2, &a[ii + a_dim1], lda);
#line 193 "zungr2.f"
	i__2 = ii + (*n - *m + ii) * a_dim1;
#line 193 "zungr2.f"
	d_cnjg(&z__2, &tau[i__]);
#line 193 "zungr2.f"
	z__1.r = 1. - z__2.r, z__1.i = 0. - z__2.i;
#line 193 "zungr2.f"
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;

/*        Set A(m-k+i,n-k+i+1:n) to zero */

#line 197 "zungr2.f"
	i__2 = *n;
#line 197 "zungr2.f"
	for (l = *n - *m + ii + 1; l <= i__2; ++l) {
#line 198 "zungr2.f"
	    i__3 = ii + l * a_dim1;
#line 198 "zungr2.f"
	    a[i__3].r = 0., a[i__3].i = 0.;
#line 199 "zungr2.f"
/* L30: */
#line 199 "zungr2.f"
	}
#line 200 "zungr2.f"
/* L40: */
#line 200 "zungr2.f"
    }
#line 201 "zungr2.f"
    return 0;

/*     End of ZUNGR2 */

} /* zungr2_ */

